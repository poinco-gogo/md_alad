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

	Output out(&fo, &sys, &ene, &atomVector);

	Integrator job(opt, &ene, &out, &atomVector);
	job.set_ptr_engine(&engine);

	const int natom = atomVector.size();
	const int nwat  = natom / 3;
	const int nfree = natom * 3 - nwat * 3;
	sys.nfree = nfree;

	const int nstep             = sys.nstep;

	cout << "REMARK Number of water molecules " << nwat << '\n';
	cout << "REMARK Degrees of freedom " << nfree << '\n';
	cout << "REMARK dt[fs] " << opt.dt_fs << '\n';
	cout << "REMARK gamma[ps-1] " << opt.langevinDamping_ps << '\n';
	cout << "REMARK T[K] " << opt.langevinTemp << '\n';

	/* generate initial velocities */
	job.reassign_velocities();

	ene.zero_force();
	ene.calc_force();
	ene.calc_kinetic_energy();
	out.print_energy(0);

	job.do_md_loop(nstep);
}
