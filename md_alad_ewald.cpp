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
#include "Lattice.hpp"
using namespace std;
class Energy
{
	public:
	double es, lj;

	Energy()
	{
		es = 0.;
		lj = 0.;
	}
};
void print_ene(int istep, Energy& ene, double K, int nfree);
bool shake(vector<Atom>& atomVector, vector<int>& shake_list);
void calc_frc(vector<Atom>& atomVector, vector<int>& lj_pair_list, vector<int>& el_pair_list, System& sys, vector<Eigen::Vector3d>& g);
double calc_kin(vector<Atom>& atomVector);
void calc_pot(vector<Atom>& atomVector, vector<int>& lj_pair_list, vector<int>& el_pair_list, System& sys, Energy& ene, vector<Eigen::Vector3d>& g);
void output(ofstream& fo, vector<Atom>& atomVector, System& sys);
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

	System sys(opt);

	mt19937 engine(static_cast<unsigned int>(sys.iseed));
	normal_distribution<double> dist(0., 1.);

	vector<Atom> atomVector;
	PSF PSFFile(argv[2], &atomVector);
	PDB PDBFile(argv[3]);
	if (!PDBFile.LoadCoords(atomVector))
		return 1;

	Energy ene;

	const int natom = atomVector.size();
	const int nwat  = natom / 3;
	const int nfree = natom * 3 - nwat * 3;

	const int print_energy_step = sys.outputEnergies;
	const int print_trj_step    = sys.DCDFreq;
	const int nstep             = sys.nstep;

        const double dt_ps = sys.dt_fs * 0.001; // in ps
        const double gamma_ps = sys.langevinDamping_ps; // in ps-1
        const double dt = dt_ps * PS2ASU;
        const double gamma = gamma_ps / PS2ASU;
	const double dt_div2 = dt * 0.5;
	const double dtdt_div2 = dt * dt * 0.5;
	const double  T = sys.langevinTemp; // K
        const double A = 1. - gamma * dt * 0.5;
        const double B = 1. + gamma * dt * 0.5;
        const double inB = 1. / B;
	
	cout << "REMARK Number of atoms " << natom << '\n';
	for (int i = 0; i < natom; i++)
	{	
		Atom& at = atomVector[i];
		at.R = sqrt(2. * BOLTZMAN * T * gamma * at.mass / dt);
	}

	cout << "REMARK Number of water molecules " << nwat << '\n';
	cout << "REMARK Degrees of freedom " << nfree << '\n';
	cout << "REMARK dt[fs] " << dt_ps * 1e3 << '\n';
	cout << "REMARK gamma[ps-1] " << gamma_ps << '\n';
	cout << "REMARK T[K] " << T << '\n';
	//cout << "REMARK Cutoff[ang.] " << system.cutoff << '\n';

	// make reciprocal vectors
	Eigen::Vector3d g1 = sys.lattice._g1();
	Eigen::Vector3d g2 = sys.lattice._g2();
	Eigen::Vector3d g3 = sys.lattice._g3();
	vector<Eigen::Vector3d> g;
	int kmax = sys.ewald_kmax;
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
	ew_self *= -sys.ewcoeff / SQRTPI * COULOMB;


	int icnt = 0;
	for (int i = 0; i < natom; i++)
	{
		Atom& at = atomVector[i];
		at.velocity.x() = dist(engine) * sqrt(BOLTZMAN * T / at.mass);
		at.velocity.y() = dist(engine) * sqrt(BOLTZMAN * T / at.mass);
		at.velocity.z() = dist(engine) * sqrt(BOLTZMAN * T / at.mass);
	}
	
	output(fo, atomVector, sys);

	vector<int> lj_pair_list, el_pair_list, shake_list;
	make_lj_pair(atomVector, lj_pair_list);
	make_el_pair(atomVector, el_pair_list);
	make_shake_pair(atomVector, shake_list);

	calc_pot(atomVector, lj_pair_list, el_pair_list, sys, ene, g);
	double Ktmp = calc_kin(atomVector);
	ene.es += ew_self;
	print_ene(0, ene, Ktmp, nfree);

	calc_frc(atomVector, lj_pair_list, el_pair_list, sys, g);
	for (int i = 0; i < atomVector.size(); i++)
	{
		Atom& at = atomVector[i];
		at.fold = at.force;
		at.rold = at.position;
		at.position = at.rold + dt * at.velocity + dtdt_div2 / at.mass * at.fold;
	}

	double accum = 0;

	for (int istep = 1; istep <= nstep; istep++)
	{
		calc_frc(atomVector, lj_pair_list, el_pair_list, sys, g);
		for (int i = 0; i < atomVector.size(); i++)
		{
			Eigen::Vector3d
			noise( dist(engine), dist(engine), dist(engine));
			Atom& at = atomVector[i];
			at.rnew = 2. * at.position - at.rold + gamma * dt_div2 * at.rold + dt * dt / at.mass * (at.force + at.R * noise);
			at.rnew = at.rnew * inB;
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
			at.vnew = 0.5 / dt * (at.rnew - at.rold);
		}

		double K = calc_kin(atomVector);
		calc_pot(atomVector, lj_pair_list, el_pair_list, sys, ene, g);
		ene.es += ew_self;
		if (istep % print_energy_step== 0)
			print_ene(istep, ene, K, nfree);

		if (istep % print_trj_step== 0)
		{
			output(fo, atomVector, sys);
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


void print_ene(int istep, Energy& ene, double K, int nfree)
{
	double totpot = ene.es + ene.lj;
	cout 
		<< setprecision(4) << fixed
		<< setw(12) << istep
		<< setw(16) << ene.lj
		<< setw(16) << ene.es
		<< setw(16) << totpot
		<< setw(16) << K
		<< setw(16) << totpot + K
		<< setw(16) << K * 2. / nfree * INVBOLTZMAN
		<< '\n';
}

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

double calc_kin(vector<Atom>& atomVector)
{
	double k = 0;
	for (int i = 0; i < atomVector.size(); i++)
	{
		Atom& at = atomVector[i];
		k += at.mass * (at.vnew).squaredNorm();
	}
	return k * 0.5;
}

void calc_frc(vector<Atom>& atomVector, vector<int>& lj_pair_list, vector<int>& el_pair_list, System& sys, vector<Eigen::Vector3d>& g)
{
	static double cutoff   = sys.cutoff;
	static double cutoff2  = cutoff * cutoff;
	static double ewcoeff  = sys.ewcoeff;
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
	static double cutoff = sys.cutoff;
	static double cutoff2 = cutoff * cutoff;
	static double boxsize = sys.box_size_x;
	static double ewcoeff = sys.ewcoeff;
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

void output(ofstream& fo, vector<Atom>& atomVector, System& sys)
{
	static const double boxsize = sys.box_size_x;
	fo << atomVector.size() << '\n' << '\n';
	fo << setprecision(6) << fixed;
	for (int i = 0; i < atomVector.size() / 3; i++) 
	{
		int j = 3 * i;
		Atom& at1 = atomVector[j];
		Atom& at2 = atomVector[j + 1];
		Atom& at3 = atomVector[j + 2];
		Eigen::Vector3d p1 = at1.position;
		Eigen::Vector3d p2 = at2.position;
		Eigen::Vector3d p3 = at3.position;
		Eigen::Vector3d vcom = (p1 + p2 + p3) / 3.;

		double d = floor(0.5 + vcom.x() / boxsize);
		p1.x() -= boxsize * d;
		p2.x() -= boxsize * d;
		p3.x() -= boxsize * d;
		
		d = floor(0.5 + vcom.y() / boxsize);
		p1.y() -= boxsize * d;
		p2.y() -= boxsize * d;
		p3.y() -= boxsize * d;
		
		d = floor(0.5 + vcom.z() / boxsize);
		p1.z() -= boxsize * d;
		p2.z() -= boxsize * d;
		p3.z() -= boxsize * d;

		fo << " " << at1.PDBAtomName[0] << setw(20) << p1.x() << setw(20) << p1.y() << setw(20) << p1.z() << '\n';
		fo << " " << at2.PDBAtomName[0] << setw(20) << p2.x() << setw(20) << p2.y() << setw(20) << p2.z() << '\n';
		fo << " " << at3.PDBAtomName[0] << setw(20) << p3.x() << setw(20) << p3.y() << setw(20) << p3.z() << '\n';
	}
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
