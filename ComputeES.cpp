#include <iostream>
#include <complex>
#include "ComputeES.hpp"
using namespace std;

ComputeES::ComputeES(const Option& opt, System& sys, vector<Atom>& atomVector, PSF& psf)
{
	this->ptr_psf            = &psf;
	this->ptr_atomVector     = &atomVector;

	this->cutoff             = opt.cutoff;

	this->usePME             = opt.usePME;
	this->ewald_kmax         = opt.ewald_kmax;
	this->ewald_tolerance    = opt.ewald_tolerance;
	this->ewcoef             = calc_ewcoef();
	this->pme_grid_x         = opt.pme_grid_x;
	this->pme_grid_y         = opt.pme_grid_y;
	this->pme_grid_z         = opt.pme_grid_z;
	this->pme_spline_order   = opt.pme_spline_order;
	this->ptr_lattice        = &sys.lattice;

	this->cutoff2            = this->cutoff * this->cutoff;
	this->ewcoef2            = this->ewcoef * this->ewcoef;

	this->boundaryType       = opt.boundaryType;

	if (boundaryType == "PBC")
	{
		if (!usePME) make_reciprocal_vectors();

		show_simulation_info_ewald();
	}
}

void ComputeES::show_simulation_info_ewald()
{
	cout << "REMARK Ewald parameters:\n"
		"REMARK     Ewald coefficient: " << ewcoef << '\n' <<
		"REMARK     Ewald tolerance  : " << ewald_tolerance << '\n' <<
		"REMARK     Ewald kmax       : " << ewald_kmax
		<< '\n';
	if (usePME)
	cout << "REMARK     use PME          : " << "yes\n" <<
		"REMARK     PME grid size x  : " << pme_grid_x << '\n' <<
		"REMARK     PME grid size y  : " << pme_grid_y << '\n' <<
		"REMARK     PME grid size z  : " << pme_grid_z << '\n' <<
		"REMARK     PME spline order : " << pme_spline_order << '\n';
}

double ComputeES::calc_ewcoef()
{
	double alpha = 1;

	while ( erfcl(alpha * cutoff) / cutoff >= ewald_tolerance )
		alpha *= 2.0;

	double alpha_lo = 0;
	double alpha_hi = alpha;

	for (int i = 0; i < 100; i++)
	{
		alpha = 0.5 * ( alpha_lo + alpha_hi );

		if ( erfcl(alpha * cutoff) / cutoff >= ewald_tolerance )
			alpha_lo = alpha;
		else
			alpha_hi = alpha;
	}

	return alpha;
}

double ComputeES::compute_force()
{
	double sum_energy = 0;

	if (this->boundaryType == "PBC")
	{
		sum_energy = compute_ewald_force();

		return sum_energy;
	}

	// NOBC case
	for (int i = 0; i < ptr_atomVector->size(); i++)
	{
		Atom& iat = ptr_atomVector->at(i);

		for (int j = i + 1; j < ptr_atomVector->size(); j++)
		{
			Atom& jat = ptr_atomVector->at(j);
			
			if (jat.checkExclusionPair(iat))
				continue;

			Eigen::Vector3d r12 = iat.position - jat.position;

			double r  = r12.norm();

			if (r > this->cutoff) continue;

			double invr = 1. / r;

			double qiqj = iat.charge * jat.charge * COULOMB;

			double e = qiqj * invr;

			sum_energy += e;

			Eigen::Vector3d f = e * invr * invr * r12;

			iat.force += f;
			jat.force -= f;
		}
	}

	return sum_energy;
}

double ComputeES::compute_ewald_force()
{
	double sum_energy = 0;

	double Udirect = 0;
	double Uself   = 0;
	double Uintra  = 0;
	double Urec    = 0;

	Lattice lattice = *ptr_lattice;

	for (int i = 0; i < 6; i++)
		tensor[i] = 0.;

	// Ewald real space summation.
	for (int i = 0; i < ptr_atomVector->size(); i++)
	{
		Atom& iat = ptr_atomVector->at(i);

		double qi    = iat.charge;
		Eigen::Vector3d ipos = iat.position;

		for (int j = i + 1; j < ptr_atomVector->size(); j++)
		{
			Atom& jat = ptr_atomVector->at(j);

			// skip if they are 1-2, 1-3 pair.
			if (iat.checkExclusionPair(jat)) continue;

			// nearest image convention.
			Eigen::Vector3d del = lattice.delta(ipos, jat.position);
			double r2 = del.squaredNorm();

			if (r2 > cutoff2) continue;
			
			double r = sqrt(r2);

			double qiqj = qi * jat.charge * COULOMB;

			double fac1 = erfcl( ewcoef * r ) / r;

			double fac2 = 2. * ewcoef / SQRTPI * exp(-ewcoef2*r2);

			Udirect  += qiqj * fac1;

			Eigen::Vector3d f = qiqj * (fac1 + fac2) / r2 * del;

			iat.force += f;
			jat.force -= f;

//			tensor[0] += f.x() * del.x();
//			tensor[2] += f.y() * del.y();
//			tensor[5] += f.z() * del.z();
		}
	}

	sum_energy += Udirect;


	// Ewald self term
	for (int i = 0; i < ptr_atomVector->size(); i++)
	{
		double qiqi = ptr_atomVector->at(i).charge;

		qiqi *= qiqi;

		Uself += qiqi;
	}
	Uself *= -ewcoef / SQRTPI * COULOMB;

	sum_energy += Uself;


	// Ewald reciprocal intra term
	for (auto& bond: ptr_psf->bondVector)
	{
		// skip dummy bond ... e.g., HT - HT
		if (bond.Kb == 0) continue;

		Atom* at1 = bond.ptr_atom1;
		Atom* at2 = bond.ptr_atom2;

		Eigen::Vector3d del = lattice.delta(at1->position, at2->position);

		double r2   = del.squaredNorm();
		double  r   = sqrt(r2);

		double fac1 = erfl(ewcoef * r) / r;
		double fac2 = 2. * ewcoef / SQRTPI * exp(-ewcoef2 * r2);

		double q1q2 = at1->charge * at2->charge * COULOMB;

		Eigen::Vector3d f   = q1q2 * (fac1 - fac2) / r2 * del;
		
		Uintra     += fac1 * q1q2;

		// intra force must be subtracted.
		at1->force -= f;
		at2->force += f;

//		tensor[0] -= f.x() * del.x();
//		tensor[2] -= f.y() * del.y();
//		tensor[5] -= f.z() * del.z();
	}

	for (auto& angle: ptr_psf->angleVector)
	{
		Atom* at1 = angle.ptr_atom1;
		Atom* at3 = angle.ptr_atom3;

		Eigen::Vector3d del = lattice.delta(at1->position, at3->position);

		double r2   = del.squaredNorm();
		double  r   = sqrt(r2);

		double fac1 = erfl(ewcoef * r) / r;
		double fac2 = 2. * ewcoef / SQRTPI * exp(-ewcoef2 * r2);

		double q1q3 = at1->charge * at3->charge * COULOMB;

		Eigen::Vector3d f   = q1q3 * (fac1 - fac2) / r2 * del;
		
		Uintra     += fac1 * q1q3;

		// intra force must be subtracted.
		at1->force -= f;
		at3->force += f;

//		tensor[0] -= f.x() * del.x();
//		tensor[2] -= f.y() * del.y();
//		tensor[5] -= f.z() * del.z();
	}

	// intra energy must be subtracted.
	sum_energy -= Uintra;

	// Ewald reciprocal sum
	if (usePME)
		sum_energy += calc_ewald_recip_pme();
	else
		sum_energy += calc_ewald_recip_direct();

	//cout << setprecision(4) << fixed;
	//cout << "DEBUG: es: " << tensor[0] << " " << tensor[2] << " " << tensor[5] << '\n';

	return sum_energy;
}

double ComputeES::calc_ewald_recip_pme()
{
	double Urec = 0;

	static helpme::Matrix<double> coordsD( ptr_atomVector->size(), 3 );
	static helpme::Matrix<double> chargeD( ptr_atomVector->size(), 1 );
	static helpme::Matrix<double> forcesD( ptr_atomVector->size(), 3 );
	static helpme::Matrix<double> virialD(1, 6);

	for (int i = 0; i < ptr_atomVector->size(); i++)
	{
		Atom& atom = ptr_atomVector->at(i);
		coordsD[i][0] = atom.position.x();
		coordsD[i][1] = atom.position.y();
		coordsD[i][2] = atom.position.z();
		chargeD[i][0] = atom.charge;
		forcesD[i][0] = 0.;
		forcesD[i][1] = 0.;
		forcesD[i][2] = 0.;
	}

	for (int i = 0; i < 6; i++)
		virialD[0][i] = 0.;

	static auto pmeD = std::unique_ptr<PMEInstanceD>(new PMEInstanceD);
	pmeD->setup(1,
			this->ewcoef,
			this->pme_spline_order,
			this->pme_grid_x,
			this->pme_grid_y,
			this->pme_grid_z,
			COULOMB,
			1
		   );
	pmeD->setLatticeVectors(
			this->ptr_lattice->_x(),
			this->ptr_lattice->_y(),
			this->ptr_lattice->_z(), 90, 90, 90,
			PMEInstanceD::LatticeType::XAligned
			);

	Urec = pmeD->computeEFVRec(0, chargeD, coordsD, forcesD, virialD);

	for (int i = 0; i < ptr_atomVector->size(); i++)
	{
		Atom& atom = ptr_atomVector->at(i);
		atom.force.x() += forcesD[i][0];
		atom.force.y() += forcesD[i][1];
		atom.force.z() += forcesD[i][2];
	}

	//cout << setprecision(4) << fixed;
	//cout << "DEBUG: pme: " << virialD[0][0] << " " << virialD[0][2] << " " << virialD[0][5] << '\n';

	this->tensor[0] += virialD[0][0];
	this->tensor[2] += virialD[0][2];
	this->tensor[5] += virialD[0][5];

	return Urec;
}

double ComputeES::calc_ewald_recip_direct()
{
	double Urec = 0;

	static complex<double> II(0., 1.);

	// loop over reciprocal lattice vectors
	for (int i = 0; i < g.size(); i++)
	{
		double g2 = g[i].squaredNorm();
		
		double dtmp;
		g[i].z() ? dtmp = 1. : dtmp = 0.5;

		double fac1 = exp(-g2 * 0.25 / ewcoef2) / g2;

		double re = 0;
		double im = 0;

		complex<double> sg(0., 0.);

		for (int j = 0; j < ptr_atomVector->size(); j++)
		{
			Atom& jat = ptr_atomVector->at(j);

			double dot = g[i].dot( jat.position );

			sg += jat.charge * exp(II * dot);

//			re += jat.charge * cos(dot);
//			im += jat.charge * sin(dot);
		}

		complex<double> conj_sg = conj(sg);

		for (int j = 0; j < ptr_atomVector->size(); j++)
		{
			Atom& jat = ptr_atomVector->at(j);

			double dot = g[i].dot( jat.position );

			complex<double> pre = II * jat.charge * exp(II * dot);

			vector< complex<double> > dsgdrj(3);

			dsgdrj[0] = pre * g[i].x();
			dsgdrj[1] = pre * g[i].y();
			dsgdrj[2] = pre * g[i].z();

			dsgdrj[0] *= conj_sg;
			dsgdrj[1] *= conj_sg;
			dsgdrj[2] *= conj_sg;

			Eigen::Vector3d f;
			double pre2 = -2. * fac1 * dtmp * 4. * PI / ptr_lattice->volume() * COULOMB;
			f.x() = pre2 * real(dsgdrj[0]);
			f.y() = pre2 * real(dsgdrj[1]);
			f.z() = pre2 * real(dsgdrj[2]);

			jat.force += f;
		}


		Urec += dtmp * fac1 * real(sg * conj(sg));
		//Urec += dtmp * fac1 * sg.squaredNorm();
	}

	Urec *= 4. * PI / ptr_lattice->volume() * COULOMB;

	return Urec;
}

double ComputeES::compute_pair_energy(Atom& at1, Atom& at2)
{
	double r = (at1.position - at2.position).norm();
	return at1.charge * at2.charge * COULOMB / r;
}

Eigen::Vector3d ComputeES::compute_pair_force(Atom& at1, Atom& at2)
{
	Eigen::Vector3d r12 = at1.position - at2.position;
	double invr = 1. / r12.norm();
	return
	at1.charge * at2.charge * COULOMB * invr * invr * invr * r12;
}

void ComputeES::make_reciprocal_vectors()
{
	Lattice lattice = *ptr_lattice;
	Eigen::Vector3d g1 = lattice._g1();
	Eigen::Vector3d g2 = lattice._g2();
	Eigen::Vector3d g3 = lattice._g3();

	const int kmax = ewald_kmax;
	int sqkmax = kmax - 1;
	sqkmax *= sqkmax;
	double dum = 0;
	for (int k = 0; k < kmax; k++)
	for (int i = -kmax + 1; i < kmax; i++)
        for (int j = -kmax + 1; j < kmax; j++)
        {
                if (k * k + j * j + i * i > sqkmax or (i==0 && j==0 && k==0))
                        continue;
                Eigen::Vector3d vtmp(i * g1.x(), j * g2.y(), k * g3.z());
                if (vtmp.norm() > dum)
                        dum = vtmp.norm();
                g.push_back(vtmp);
        }
        cout << "REMARK gmax = " << dum << '\n';
}
