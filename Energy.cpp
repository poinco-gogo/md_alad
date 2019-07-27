#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include "Energy.hpp"
#include "common.hpp"
using namespace std;

Energy::Energy(const Option& opt, System& sys, vector<Atom>& atomVector, PSF& psf)
{
	this->ptr_sys            = &sys;
	this->ptr_atomVector     = &atomVector;

	this->outputEnergies     = opt.outputEnergies;

	this->cutoff             = opt.cutoff;

	this->usePME             = opt.usePME;
	this->ewald_kmax         = opt.ewald_kmax;
	this->pme_grid_x         = opt.pme_grid_x;
	this->pme_grid_y         = opt.pme_grid_y;
	this->pme_grid_z         = opt.pme_grid_z;
	this->pme_spline_order   = opt.pme_spline_order;
	this->ewald_tolerance    = opt.ewald_tolerance;

	ComputeLJ tmp_lj(opt, sys, atomVector);
	ComputeES tmp_es(opt, sys, atomVector, psf);

	this->vlj = tmp_lj;
	this->ves = tmp_es;

	if (opt.boundaryType == "PBC") make_reciprocal_vectors();

	make_lj_pair_list();
	make_el_pair_list();
}

void Energy::show_simulation_info()
{
	cout 
		<< "REMARK CUTOFF                 " << this->cutoff
		<< '\n';
}

void Energy::calc_kinetic_energy()
{
	this->kinetic = 0;

	for (auto& atom: *ptr_atomVector)
	{
		this->kinetic += atom.mass * atom.vnew.squaredNorm();
	}

	this->kinetic *= 0.5;
}

void Energy::make_lj_pair_list()
{
	int icnt = 0;
	for (auto& at: *ptr_atomVector)
	{
		if (at.PDBAtomName == "OH2")
			lj_pair_list.push_back(icnt);
		++icnt;
	}
}

void Energy::make_el_pair_list()
{
	for (int i = 0; i < ptr_atomVector->size(); i++)
	{
		Atom& at1 = ptr_atomVector->at(i);

		for (int j = i + 1; j < ptr_atomVector->size(); j++)
		{
			Atom& at2 = ptr_atomVector->at(j);

			if (i / 3 == j / 3) continue;

			el_pair_list.push_back(i);
			el_pair_list.push_back(j);
		}
	}
}

void Energy::make_reciprocal_vectors()
{
	Eigen::Vector3d g1 = ptr_sys->lattice._g1();
	Eigen::Vector3d g2 = ptr_sys->lattice._g2();
	Eigen::Vector3d g3 = ptr_sys->lattice._g3();

	int kmax = ewald_kmax;

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

void Energy::calc_force()
{
	lj = vlj.compute_force();
	es = ves.compute_force();
/*	Lattice& lattice = ptr_sys->lattice;
	vector<Atom>& atomVector = *ptr_atomVector;

	//static double cutoff   = ene.cutoff;
	static double cutoff2  = cutoff * cutoff;
	//static double ewcoeff  = ene.ewcoeff;
	static double ewcoeff2 = ewcoeff * ewcoeff;
	static double boxsize  = lattice._x();
	static double factor = 1. / (4. * ewcoeff * ewcoeff);
	static double const_intra   = -2. * ewcoeff / SQRTPI;
	static double const_recipro = 4. * PI / lattice.volume();
	static double A = 582. * 1e3;
	static double A_2 = A * 2.;
	static double B = 595.0;
	for (int i = 0; i < atomVector.size(); i++)
	{
		atomVector[i].force = V3ZERO;
	}

	for (int i = 0; i < lj_pair_list.size(); i++)
	{
		Atom& at1 = atomVector[lj_pair_list[i]];

		for (int j = i + 1; j < lj_pair_list.size(); j++)
		{
			Atom& at2 = atomVector[lj_pair_list[j]];
			Eigen::Vector3d del
			= lattice.delta(at1.position ,at2.position);
			double r2 = (del).squaredNorm();
			if (r2 > cutoff2) continue;
			double r6 = r2 * r2 * r2;
			double r8 = r6 * r2;
			Eigen::Vector3d f12 = 6. / r8 * (B - A_2 / r6) * del;
			at1.force -= f12;
			at2.force += f12;
		}
	//	if (!i) print(f12);
	}
//	print(atomVector[0].fnew);

	// ewald intra force
	for (int i = 0; i < atomVector.size() / 3; i++)
	{
		Eigen::Vector3d frc = V3ZERO;
		int j = 3 * i;
		Atom& at1 = atomVector[j];
		Atom& at2 = atomVector[j + 1];
		Atom& at3 = atomVector[j + 2];
		Eigen::Vector3d del1
			= lattice.delta(at1.position ,at2.position);
		double adel1 = del1.norm();
		double sqadel1 = adel1 * adel1;
		frc = at1.charge * at2.charge * COULOMB *
			(const_intra * exp(-ewcoeff2 * sqadel1)
			 + erfl(ewcoeff * adel1) / adel1) / sqadel1 * del1;
		at1.force -= frc;
		at2.force += frc;
		Eigen::Vector3d del2
			= lattice.delta(at1.position ,at3.position);
		double adel2 = del2.norm();
		double sqadel2 = adel2 * adel2;
		frc = at1.charge * at3.charge * COULOMB *
			(const_intra * exp(-ewcoeff2 * sqadel2)
			 + erfl(ewcoeff * adel2) / adel2) / sqadel2 * del2;
		at1.force -= frc;
		at3.force += frc;
		Eigen::Vector3d del3
			= lattice.delta(at3.position ,at2.position);
		double adel3 = del3.norm();
		double sqadel3 = adel3 * adel3;
		frc = at3.charge * at2.charge * COULOMB *
			(const_intra * exp(-ewcoeff2 * sqadel3)
			 + erfl(ewcoeff * adel3) / adel3) / sqadel3 * del3;
		at3.force -= frc;
		at2.force += frc;
	}


	// ewald direct space force summation
	for (int i = 0; i < el_pair_list.size() / 2; i++)
	{
		Atom& at1 = atomVector[el_pair_list[2 * i]];
		Atom& at2 = atomVector[el_pair_list[2 * i + 1]];

		Eigen::Vector3d del
			= lattice.delta(at1.position ,at2.position);
		double r2 = del.squaredNorm();
		if (r2 > cutoff2) continue;
		double r = sqrt(r2);
		Eigen::Vector3d frc = V3ZERO;
		frc = at1.charge * at2.charge * COULOMB *
			(-const_intra * exp(-ewcoeff2 * r2)
			+ erfc(ewcoeff*r) / r) / r2 *del;
		at1.force += frc;
		at2.force -= frc;
	}


	// ewald reciprocal space summation
	for (int ii = 0; ii < atomVector.size(); ii++)
	{
		Atom& iat = atomVector[ii];
		Eigen::Vector3d frc = V3ZERO;
	for (int i = 0; i < g.size(); i++)
	{
		double ag = g[i].norm();
		double agag = ag * ag;
		double pre, dtmp;
		pre = 0;
		for (int j = 0; j < atomVector.size(); j++)
		{
			Atom& jat = atomVector[j];
			Eigen::Vector3d del = iat.position - jat.position;
			del.x() -= boxsize * floor(del.x() / boxsize + 0.5);
			del.y() -= boxsize * floor(del.y() / boxsize + 0.5);
			del.z() -= boxsize * floor(del.z() / boxsize + 0.5);
			double dot = g[i].dot(del);

			pre += jat.charge * sin(dot);
		}
		g[i].z() ? dtmp = 1. : dtmp = 0.5;
		frc = frc + (dtmp * exp(-agag * factor) / agag * pre) * g[i];
	}
		iat.force += frc * const_recipro * COULOMB * iat.charge;
//		cout << vabs(frc) << '\n';
	}
//	print(atomVector[0].fnew);*/
}

void Energy::calc_potential_energy()
{
	Lattice& lattice = ptr_sys->lattice;
	vector<Atom>& atomVector = *ptr_atomVector;

	//static double cutoff = ene.cutoff;
	static double cutoff2 = cutoff * cutoff;
	static double boxsize = lattice._x();
	//static double ewcoeff = ene.ewcoeff;
	static double factor = 1. / (4. * ewcoeff * ewcoeff);
	static double const_recipro = 4. * PI / lattice.volume();
	static double A = 582. * 1e3;
	static double B = 595.0;
	static double qO = -0.834;
	static double qH = 0.417;

	lj = 0;
	es = 0;
	double ew_direct  = 0;
	double ew_recipro = 0;
	double ew_intra   = 0;
	double U = 0;

	for (int i = 0; i < lj_pair_list.size(); i++)
	{
		Atom& iat = ptr_atomVector->at( lj_pair_list[i] );
		for (int j = i + 1; j < lj_pair_list.size(); j++)
		{
			Atom& jat = ptr_atomVector->at( lj_pair_list[j] );
			Eigen::Vector3d del
			= lattice.delta(iat.position ,jat.position);
			double roo2 = del.squaredNorm();
			if (roo2 > cutoff2) continue;
			double roo6 = roo2 * roo2 * roo2;
			lj += A / (roo6 * roo6) - B / roo6;
		}
	}

	// ewald intra energy
	for (int i = 0; i < atomVector.size() / 3; i++)
	{
		int j = 3 * i;
		Atom& at1 = atomVector[j];
		Atom& at2 = atomVector[j + 1];
		Atom& at3 = atomVector[j + 2];
		Eigen::Vector3d del1
			= lattice.delta(at1.position ,at2.position);
		double adel1 = del1.norm();
		Eigen::Vector3d del2
			= lattice.delta(at1.position ,at3.position);
		double adel2 = del2.norm();
		Eigen::Vector3d del3
			= lattice.delta(at3.position ,at2.position);
		double adel3 = del3.norm();
		ew_intra += at1.charge*at2.charge*erfl(ewcoeff*adel1)/adel1;
		ew_intra += at1.charge*at3.charge*erfl(ewcoeff*adel2)/adel2;
		ew_intra += at3.charge*at2.charge*erfl(ewcoeff*adel3)/adel3;
	}

	// ewald direct space summation
	for (int i = 0; i < el_pair_list.size() / 2; i++)
	{
		Atom& at1 = atomVector[el_pair_list[2 * i]];
		Atom& at2 = atomVector[el_pair_list[2 * i + 1]];

		Eigen::Vector3d del
			= lattice.delta(at1.position ,at2.position);
		double r = del.norm();
		if (r > cutoff) continue;
		ew_direct += at1.charge * at2.charge * erfc(ewcoeff*r)/r;
	}

	// ewald reciprocal space summation
	for (int i = 0; i < g.size(); i++)
	{
		double ag = g[i].norm();
		double agag = ag * ag;
		double re, im, dtmp;
		re = 0;
		im = 0;
		for (int j = 0; j < atomVector.size(); j++)
		{
			Eigen::Vector3d del = atomVector[j].position;
			del.x() -= boxsize * floor(del.x() / boxsize + 0.5);
			del.y() -= boxsize * floor(del.y() / boxsize + 0.5);
			del.z() -= boxsize * floor(del.z() / boxsize + 0.5);
			double dot = g[i].dot(del);

			re += atomVector[j].charge * cos(dot);
			im += atomVector[j].charge * sin(dot);
		}
		g[i].z() ? dtmp = 1. : dtmp = 0.5;
		ew_recipro += dtmp * exp(-agag * factor) / agag * (re * re + im * im);
	}

	ew_recipro *= const_recipro;

//	cout << "el_direct  " << el_direct * 332.0636 << '\n';
//	cout << "el_intra    " << el_intra * 332.0636 << '\n';
//	cout << "el_recipro " << el_recipro * 332.0636<< '\n';

	es += ew_direct + ew_recipro - ew_intra;
	es *= COULOMB;
}

void Energy::zero_force()
{
	for (auto& atom: *ptr_atomVector)
		atom.force = V3ZERO;
}
