#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <random>
#include "common.hpp"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Atom.hpp"
#include "PSF.hpp"
#include "PDB.hpp"
#include "LoadParm.hpp"
#include "Option.hpp"
#include "System.hpp"
#include "Energy.hpp"
#include "Integrator.hpp"
#include "Output.hpp"
#include "Lattice.hpp"
using namespace std;
bool shake(vector<Atom>& atomVector, vector<int>& shake_list);
void make_shake_pair(vector<Atom>& atomVector, vector<int>& shake_list);
int main (int argc, char** argv)
{
	if (argc < 3)
	{
		cout << "usage: ./a.out trjout config\n";
		return 1;
	}
	
	ofstream fo(argv[1]);
	Option opt(argv[2]);
	if (!opt.load_config())
		return 1;

	System     sys(opt);

	mt19937 engine(static_cast<unsigned int>(sys.iseed));
	normal_distribution<double> dist(0., 1.);

	vector<Atom> atomVector;
	PSF PSFFile(sys.structure, &atomVector);
	PDB PDBFile(sys.coordinates);
	if (!PDBFile.LoadCoords(atomVector))
		return 1;

	LoadParm ALL22(sys.parameters);
	if (!PSFFile.set_bond_parm(ALL22.bondParmVector)) return 0;
	if (!PSFFile.set_angle_parm(ALL22.angleParmVector)) return 0;
	if (!PSFFile.set_lj_parm(ALL22.LJParmVector)) return 0;

	PSFFile.make_exclusion_vector();

	Energy     ene(opt, sys, atomVector, PSFFile);
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

	double ew_self = 0;
	for (int i = 0; i < atomVector.size(); i++)
	{	
		ew_self += atomVector[i].charge * atomVector[i].charge;
	}
	ew_self *= -ene.ewcoeff / SQRTPI * COULOMB;

	/* generate initial velocities */
	job.reassign_velocities();

	vector<int> shake_list;
	make_shake_pair(atomVector, shake_list);

	ene.zero_force();
	ene.calc_force();
	ene.calc_kinetic_energy();
	out.print_energy(0);

	job.initial_posi_velret();

	for (int istep = 1; istep <= nstep; istep++)
	{
		ene.zero_force();
		ene.calc_force();
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
		//ene.calc_potential_energy();
		//ene.es += ew_self;
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
