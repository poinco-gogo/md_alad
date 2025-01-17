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

	this->rigidBonds         = opt.rigidBonds;
	this->rigidIndexes       = opt.rigidIndexes;
	this->rigidTolerance     = opt.rigidTolerance;
	this->rigidIterations    = opt.rigidIterations;

	this->boundaryType       = opt.boundaryType;
	this->wrapAll            = opt.wrapAll;
	this->box_size_x         = opt.box_size_x;
	this->box_size_y         = opt.box_size_y;
	this->box_size_z         = opt.box_size_z;

	this->iseed              = opt.iseed;

	if (this->boundaryType == "PBC")
	{
		make_lattice(box_size_x, box_size_y, box_size_z);
	}

	this->nfree = 0;

	show_simulation_info();
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

	if (this->binvelocities.size())
	{
		cout
		<< "REMARK RESTART VELOCITIES     " << this->binvelocities
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

	if (this->wrapAll)
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
