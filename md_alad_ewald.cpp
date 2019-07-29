#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <random>
#include <chrono>
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
#include "NAMDBin.hpp"
#include "Lattice.hpp"
#include "Random.hpp"
using namespace std;
int main (int argc, char** argv)
{
	if (argc < 2)
	{
		cout << "usage: ./a.out config\n";
		return 1;
	}
	
	Option opt(argv[1]);
	if (!opt.load_config())
		return 1;

	System     sys(opt);

	Random random(sys.iseed);

	vector<Atom> atomVector;
	PSF PSFFile(sys.structure, &atomVector);
	PDB PDBFile(sys.coordinates);
	if (sys.bincoordinates.size())
	{
		NAMDBin bincoor(sys.bincoordinates.c_str(), "coor");
		if (!bincoor.read_fi(atomVector)) return 1;
	}
	else if (!PDBFile.LoadCoords(atomVector)) return 1;

	if (opt.rigidBonds) PSFFile.make_water_shake_bond_array();

	LoadParm ALL22(sys.parameters);
	if (!PSFFile.set_bond_parm(ALL22.bondParmVector)) return 0;
	if (!PSFFile.set_angle_parm(ALL22.angleParmVector)) return 0;
	if (!PSFFile.set_lj_parm(ALL22.LJParmVector)) return 0;

	PSFFile.make_exclusion_vector();

	Energy     ene(opt, sys, atomVector, PSFFile);

	Output out(opt, &sys, &ene, &atomVector);

	Integrator job(opt, &ene, &out, &atomVector, &random);

	sys.nfree = 3 * atomVector.size() - ene.vbnd._num_shake_bond();

	cout << "REMARK Degrees of freedom " << sys.nfree << '\n';
	cout << "REMARK dt[fs] " << opt.dt_fs << '\n';
	cout << "REMARK gamma[ps-1] " << opt.langevinDamping_ps << '\n';
	cout << "REMARK T[K] " << opt.langevinTemp << '\n';

	/*  velocities at zero step */
	if (sys.binvelocities.size())
	{
		NAMDBin binvel(sys.binvelocities.c_str(), "vel");
		if (!binvel.read_fi(atomVector)) return 1;
	}
	else job.reassign_velocities();

	ene.zero_force();
	ene.calc_force();
	ene.calc_kinetic_energy();
	out.print_energy(0);

	auto start = std::chrono::system_clock::now();
	job.do_md_loop(sys.nstep);
	auto end   = std::chrono::system_clock::now();
	auto dur   = end - start;
	double msec
	= std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
	cout << "REMARK ----------------------------------\n"
		"REMARK elapsed:  " << msec / 1000 << " sec.\n\n";

	out.output_namdbin("coor");
	out.output_namdbin("vel");

	out.close_dcd();
}
