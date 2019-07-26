#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <ctime>
#include <random>
#include "common.hpp"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Atom.hpp"
#include "PSF.hpp"
#include "PDB.hpp"
#include "Option.hpp"
#include "System.hpp"
#include "Energy.hpp"
#include "Integrator.hpp"
#include "Output.hpp"
#include "Lattice.hpp"
using namespace std;
bool shake(vector<Atom>& atomVector, vector<int>& shake_list);
void calc_frc(vector<Atom>& atomVector, vector<int>& lj_pair_list, vector<int>& el_pair_list, System& sys, Energy& ene, vector<Eigen::Vector3d>& g);
void calc_pot(vector<Atom>& atomVector, vector<int>& lj_pair_list, vector<int>& el_pair_list, System& sys, Energy& ene, vector<Eigen::Vector3d>& g);
void make_lj_pair(vector<Atom>& atomVector, vector<int>& lj_pair_list);
void make_el_pair(vector<Atom>& atomVector, vector<int>& el_pair_list);
void make_shake_pair(vector<Atom>& atomVector, vector<int>& shake_list);
int main (int argc, char** argv)
{
	if (argc < 4)
	{
		cout << "usage: ./a.out trjout psf pdb config\n";
		return 1;
	}
	
	ofstream fo(argv[1]);
	Option opt(argv[4]);
	if (!opt.load_config())
		return 1;

	System     sys(opt);

	mt19937 engine(static_cast<unsigned int>(sys.iseed));
	normal_distribution<double> dist(0., 1.);

	vector<Atom> atomVector;
	PSF PSFFile(argv[2], &atomVector);
	PDB PDBFile(argv[3]);
	if (!PDBFile.LoadCoords(atomVector))
		return 1;

	Energy     ene(opt, &atomVector);
	Integrator job(opt, &atomVector);
	job.set_ptr_engine(&engine);

	Output out(&atomVector, &sys, &ene);

	const int natom = atomVector.size();
	const int nwat  = natom / 3;
	const int nfree = natom * 3 - nwat * 3;
	sys.nfree = nfree;

	const int print_energy_step = ene.outputEnergies;
	const int print_trj_step    = sys.DCDFreq;
	const int nstep             = sys.nstep;

	cout << "REMARK Number of atoms " << natom << '\n';
	cout << "REMARK Number of water molecules " << nwat << '\n';
	cout << "REMARK Degrees of freedom " << nfree << '\n';
	cout << "REMARK dt[fs] " << opt.dt_fs << '\n';
	cout << "REMARK gamma[ps-1] " << opt.langevinDamping_ps << '\n';
	cout << "REMARK T[K] " << opt.langevinTemp << '\n';

	// make reciprocal vectors
	Eigen::Vector3d g1 = sys.lattice._g1();
	Eigen::Vector3d g2 = sys.lattice._g2();
	Eigen::Vector3d g3 = sys.lattice._g3();
	vector<Eigen::Vector3d> g;
	int kmax = ene.ewald_kmax;
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

	double ew_self = 0;
	for (int i = 0; i < atomVector.size(); i++)
	{	
		ew_self += atomVector[i].charge * atomVector[i].charge;
	}
	ew_self *= -ene.ewcoeff / SQRTPI * COULOMB;

	/* generate initial velocities */
	job.reassign_velocities();

	vector<int> lj_pair_list, el_pair_list, shake_list;
	make_lj_pair(atomVector, lj_pair_list);
	make_el_pair(atomVector, el_pair_list);
	make_shake_pair(atomVector, shake_list);

	calc_pot(atomVector, lj_pair_list, el_pair_list, sys, ene, g);
	ene.calc_kinetic_energy();
	ene.es += ew_self;
	out.print_energy(0);

	calc_frc(atomVector, lj_pair_list, el_pair_list, sys, ene, g);

	job.initial_posi_velret();

	for (int istep = 1; istep <= nstep; istep++)
	{
		calc_frc(atomVector, lj_pair_list, el_pair_list, sys, ene, g);
		for (int i = 0; i < atomVector.size(); i++)
		{
			Eigen::Vector3d
			noise( dist(engine), dist(engine), dist(engine));
			Atom& at = atomVector[i];
			at.rnew = 2. * at.position - at.rold + job.gamma * job.dt_div2 * at.rold + job.dt * job.dt * at.invmass * (at.force + at.R * noise);
			at.rnew = at.rnew * job.inB;
		}

		for (int ishake = 0; ishake < 100; ishake++)
		{
			if (shake(atomVector, shake_list))
				break;
			if (ishake == 99)
			{
				cerr << "error: shake does not converged.\n";
				return 1;
			}
		}

		for (int i = 0; i < atomVector.size(); i++)
		{
			Atom& at = atomVector[i];
			at.vnew = 0.5 / job.dt * (at.rnew - at.rold);
		}

		ene.calc_kinetic_energy();
		calc_pot(atomVector, lj_pair_list, el_pair_list, sys, ene, g);
		ene.es += ew_self;
		if (istep % print_energy_step== 0)
			out.print_energy(istep);

		if (istep % print_trj_step== 0)
		{
			out.output_xyz(fo);
		}

		for (int i = 0; i < atomVector.size(); i++)
		{
			Atom& at = atomVector[i];
			at.rold = at.position;
			at.position = at.rnew;
			at.velocity = at.vnew;
		}
	}
}
/////////////////////////  end of main program

bool shake(vector<Atom>& atomVector, vector<int>& shake_list)
{
	static const double eps = 1e-6;
	static const double eps2 = eps * eps;
	static const double rOH = 0.9572;
	static const double aHOH = 104.52 * DEG2RAD;
	static const double dOH = rOH * rOH;
	static const double rHH = rOH * sin(aHOH/2.) * 2.;
	static const double dHH = rHH * rHH;
	for (int i = 0; i < shake_list.size() / 2; i++)
	{
		Atom& at1 = atomVector[shake_list[2 * i]];
		Atom& at2 = atomVector[shake_list[2 * i + 1]];
		double gamma;
		if (at1.PDBAtomName[0] == 'O' && at2.PDBAtomName[0] == 'H')
		{
			gamma = (dOH - (at1.rnew - at2.rnew).squaredNorm()) /
				(2.*(1./at1.mass+1./at2.mass)*((at1.position - at2.position).dot(at1.rnew - at2.rnew)));
		}
		else
		{
			gamma = (dHH - (at1.rnew - at2.rnew).squaredNorm()) /
				(2.*(1./at1.mass+1./at2.mass)*((at1.position - at2.position).dot(at1.rnew - at2.rnew)));

		}
		at1.rnew = at1.rnew + gamma * (at1.position - at2.position) / at1.mass;
		at2.rnew = at2.rnew + gamma * (at2.position - at1.position) / at2.mass;
	}

	for (int i = 0; i < shake_list.size() / 2; i++)
	{
		Atom& at1 = atomVector[shake_list[2 * i]];
		Atom& at2 = atomVector[shake_list[2 * i + 1]];
		double r = (at1.rnew - at2.rnew).norm();
		double error;
		if (at1.PDBAtomName[0]  == 'O' && at2.PDBAtomName[0] == 'H')
		{
			error = abs(r - rOH);
		}
		else
		{
			error = abs(r - rHH);
		}
		if (error > eps) return false;
	}

	return true;
}

void calc_frc(vector<Atom>& atomVector, vector<int>& lj_pair_list, vector<int>& el_pair_list, System& sys, Energy& ene, vector<Eigen::Vector3d>& g)
{
	static double cutoff   = ene.cutoff;
	static double cutoff2  = cutoff * cutoff;
	static double ewcoeff  = ene.ewcoeff;
	static double ewcoeff2 = ewcoeff * ewcoeff;
	static double boxsize  = sys.box_size_x;
	static double factor = 1. / (4. * ewcoeff * ewcoeff);
	static double const_intra   = -2. * ewcoeff / SQRTPI;
	static double const_recipro = 4. * PI / sys.lattice.volume();
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
			= sys.lattice.delta(at1.position ,at2.position);
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
			= sys.lattice.delta(at1.position ,at2.position);
		double adel1 = del1.norm();
		double sqadel1 = adel1 * adel1;
		frc = at1.charge * at2.charge * COULOMB *
			(const_intra * exp(-ewcoeff2 * sqadel1) 
			 + erfl(ewcoeff * adel1) / adel1) / sqadel1 * del1;
		at1.force -= frc;
		at2.force += frc;
		Eigen::Vector3d del2
			= sys.lattice.delta(at1.position ,at3.position);
		double adel2 = del2.norm();
		double sqadel2 = adel2 * adel2;
		frc = at1.charge * at3.charge * COULOMB *
			(const_intra * exp(-ewcoeff2 * sqadel2) 
			 + erfl(ewcoeff * adel2) / adel2) / sqadel2 * del2;
		at1.force -= frc;
		at3.force += frc;
		Eigen::Vector3d del3
			= sys.lattice.delta(at3.position ,at2.position);
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
			= sys.lattice.delta(at1.position ,at2.position);
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
//	print(atomVector[0].fnew);
}

void calc_pot(vector<Atom>& atomVector, vector<int>& lj_pair_list, vector<int>& el_pair_list, System& sys, Energy& ene, vector<Eigen::Vector3d>& g)
{
	static double cutoff = ene.cutoff;
	static double cutoff2 = cutoff * cutoff;
	static double boxsize = sys.box_size_x;
	static double ewcoeff = ene.ewcoeff;
	static double factor = 1. / (4. * ewcoeff * ewcoeff);
	static double const_recipro = 4. * PI / sys.lattice.volume();
	static double A = 582. * 1e3;
	static double B = 595.0;
	static double qO = -0.834;
	static double qH = 0.417;

	ene.lj = 0;
	ene.es = 0;
	double ew_direct  = 0;
	double ew_recipro = 0;
	double ew_intra   = 0;
	double U = 0;

	for (int i = 0; i < lj_pair_list.size(); i++)
	{
		Atom& iat = atomVector[lj_pair_list[i]];
		for (int j = i + 1; j < lj_pair_list.size(); j++)
		{
			Atom& jat = atomVector[lj_pair_list[j]];
			Eigen::Vector3d del
			= sys.lattice.delta(iat.position ,jat.position);
			double roo2 = del.squaredNorm();
			if (roo2 > cutoff2) continue;
			double roo6 = roo2 * roo2 * roo2;
			ene.lj += A / (roo6 * roo6) - B / roo6;
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
			= sys.lattice.delta(at1.position ,at2.position);
		double adel1 = del1.norm();
		Eigen::Vector3d del2
			= sys.lattice.delta(at1.position ,at3.position);
		double adel2 = del2.norm();
		Eigen::Vector3d del3
			= sys.lattice.delta(at3.position ,at2.position);
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
			= sys.lattice.delta(at1.position ,at2.position);
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

	ene.es += ew_direct + ew_recipro - ew_intra;
	ene.es *= COULOMB;
}

void make_lj_pair(vector<Atom>& atomVector, vector<int>& lj_pair_list)
{
	for (int i = 0; i < atomVector.size(); i++)
	{
		if (atomVector[i].PDBAtomName == "OH2")
			lj_pair_list.push_back(i);
	}
}

void make_el_pair(vector<Atom>& atomVector, vector<int>& el_pair_list)
{
	for (int i = 0; i < atomVector.size(); i++)
	{
		Atom& at1 = atomVector[i];

		for (int j = i + 1; j < atomVector.size(); j++)
		{
			Atom& at2 = atomVector[j];

			if (i / 3 == j / 3) continue;

			el_pair_list.push_back(i);
			el_pair_list.push_back(j);
		}
	}
}

void make_shake_pair(vector<Atom>& atomVector, vector<int>& shake_list)
{
	for (int i = 0; i < atomVector.size() / 3; i++)
	{
		int j = 3 * i;
		shake_list.push_back(j);
		shake_list.push_back(j+1);
		shake_list.push_back(j);
		shake_list.push_back(j+2);
		shake_list.push_back(j+1);
		shake_list.push_back(j+2);
	}
}
