#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include "System.hpp"
#include "common.hpp"
using namespace std;

System::System(const Option& opt)
{
	this->structure          = opt.structure;
	this->coordinates        = opt.coordinates;
	this->parameters         = opt.parameters;
	this->topparfile         = opt.topparfile;
	this->outputname         = opt.outputname;
	this->bincoordinates     = opt.bincoordinates;
	this->binvelocities      = opt.binvelocities;
	this->nstep              = opt.nstep;
	this->firsttimestep      = opt.firsttimestep;
	this->DCDFreq            = opt.DCDFreq;
	this->outputEnergies     = opt.outputEnergies;

	this->initialTemp        = opt.initialTemp;

	this->rigidBonds         = opt.rigidBonds;
	this->rigidIndexes       = opt.rigidIndexes;
	this->rigidTolerance     = opt.rigidTolerance;
	this->rigidIterations    = opt.rigidIterations;

	this->boundaryType       = opt.boundaryType;
	this->wrapAll            = opt.wrapAll;
	this->box_size_x         = opt.box_size_x;
	this->box_size_y         = opt.box_size_y;
	this->box_size_z         = opt.box_size_z;

	this->cutoff             = opt.cutoff;

	this->usePME             = opt.usePME;
	this->ewald_kmax         = opt.ewald_kmax;
	this->pme_grid_x         = opt.pme_grid_x;
	this->pme_grid_y         = opt.pme_grid_y;
	this->pme_grid_z         = opt.pme_grid_z;
	this->pme_spline_order   = opt.pme_spline_order;
	this->ewald_tolerance    = opt.ewald_tolerance;

	this->iseed              = opt.iseed;

	this->dt_fs              = opt.dt_fs;

	this->langevinTemp       = opt.langevinTemp;
	this->langevinDamping_ps = opt.langevinDamping_ps;
	this->langevin           = opt.langevin;

	if (this->boundaryType == "PBC")
	{
		make_lattice(box_size_x, box_size_y, box_size_z);
	}
}

void System::make_lattice(double x, double y, double z)
{
	Lattice ltmp(x, y, z);
	this->lattice = ltmp;
}

void System::show_simulation_info()
{
	cout
                << "REMARK SEED VALUE             " << this->iseed
                << '\n'
                << "REMARK NUMBER OF STEPS        " << this->nstep
                << '\n';

	cout 
		<< "REMARK CUTOFF                 " << this->cutoff
		<< '\n';

	cout
                << "REMARK COORDINATE PDB         " << this->coordinates
                << '\n'
                << "REMARK STRUCTURE FILE         " << this->structure
                << '\n'
                << "REMARK PARAMETERS             " << this->parameters
                << '\n';

	if (this->topparfile != "")
		cout << "REMARK TOPPARFILE             "
			<< this->topparfile << '\n';
	
	if (this->bincoordinates.size())
	{
		cout
		<< "REMARK RESTART CONFIGURATION  " << this->bincoordinates
		<< '\n'
		<< "REMARK COORDINATE PDB WILL BE IGNORED."
		<< '\n';
	}

	if (this->boundaryType == "PBC")
	{
		cout << "REMARK Periodic boundary condition:\n"
			"REMARK     box_size_x: " << this->box_size_x << '\n' <<
			"REMARK     box_size_y: " << this->box_size_y << '\n' <<
			"REMARK     box_size_z: " << this->box_size_z << '\n';
	}
	else
	{
		cout << "REMARK NO Periodic boundary condition\n";
	}

	if (this->wrapAll == "on")
	{
		cout << "REMARK WRAPPING ALL CLUSTERS AROUND PERIODIC BOUNDARIES ON OUTPUT.\n";
	}
}

void System::make_solute_index(vector<Atom>& atomVector)
{
	for (auto& atom: atomVector)
	{
		if (atom.PSFResName != "TIP3")
			soluteIndex.push_back(atom.PSFIndex);
	}
}

void System::make_water_index(vector<Atom>& atomVector)
{
	for (auto& atom: atomVector)
	{
		if (atom.PDBAtomName == "OH2")
			waterIndex.push_back(atom.PSFIndex);
	}
}
