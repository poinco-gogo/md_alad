#include <boost/algorithm/string.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <random>
#include <cmath>
#include "System.hpp"
#include "common.hpp"
using namespace std;

System::System(std::string filename)
{
	this->filename = filename;
	reset();
}

void System::reset()
{
	firsttimestep = 0;
	com_remove = "off";
	cutoff = 9999.;
        initialTemp = -1;
        nstep = 0;
	rigidBonds      = "none";
        rigidTolerance  = 1e-8;
        rigidIterations = 100;
	rigidIndexes    = "";
	iseed          = 1234;
	boundaryType   = "NOBC";
	box_size_x     = -1;
	box_size_y     = -1;
	box_size_z     = -1;
	ewald_kmax     = 20;
	ewald_tolerance = 1e-6;
	usePME         = "no";
	pme_grid_x     = 32;
	pme_grid_y     = 32;
	pme_grid_z     = 32;
	pme_spline_order = 4;
	wrapAll        = "off";
}

bool System::load_config()
{
	ifstream fc(filename.c_str());

	std::string s;
	while (getline(fc, s))
        {
                istringstream is(s);
		std::string stmp;
                if (s.find("#", 0) != string::npos || s.empty())
                {
                        // found comment line.
                        continue;
                }
		else if (s.find("iseed", 0) != string::npos)
                {
                        is >> stmp >> this->iseed;
                }
                else if (s.find("temperature", 0) != string::npos)
                {
                        is >> stmp >> this->initialTemp;
                }
                else if (s.find("numstep", 0) != string::npos)
                {
                        is >> stmp >> this->nstep;
                }
                else if (s.find("firsttimestep", 0) != string::npos)
                {
                        is >> stmp >> this->firsttimestep;
                }
                else if (s.find("DCDFreq", 0) != string::npos)
                {
                        is >> stmp >> this->DCDFreq;
                }
		else if (s.find("outputEnergies", 0) != string::npos)
                {
                        is >> stmp >> this->outputEnergies;
                }
                else if (s.find("structure", 0) != string::npos)
                {
                        is >> stmp >> this->structure;
                }
                else if (s.find("bincoordinates", 0) != string::npos)
                {
                        is >> stmp >> this->bincoordinates;
                }
                else if (s.find("binvelocities", 0) != string::npos)
                {
                        is >> stmp >> this->binvelocities;
                }
                else if (s.find("coordinates", 0) != string::npos)
                {
                        is >> stmp >> this->coordinates;
                }
                else if (s.find("restartfilename", 0) != string::npos)
                {
                        is >> stmp >> this->restartfilename;
                }
                else if (s.find("parameters", 0) != string::npos)
                {
                        is >> stmp >> this->parameters;
		}
                else if (s.find("topparfile", 0) != string::npos)
                {
                        is >> stmp >> this->topparfile;
		}
                else if (s.find("outputname", 0) != string::npos)
                {
                        is >> stmp >> this->outputname;
		}
                else if (s.find("rigidBonds", 0) != string::npos)
                {
                        is >> stmp >> this->rigidBonds;
                }
		else if (s.find("rigidTolerance", 0) != string::npos)
                {
                        is >> stmp >> this->rigidTolerance;
                }
                else if (s.find("rigidIterations", 0) != string::npos)
                {
                        is >> stmp >> this->rigidIterations;
                }
                else if (s.find("rigidIndexes", 0) != string::npos)
                {
                        is >> stmp >> this->rigidIndexes;
                }
                else if (s.find("com_remove", 0) != string::npos)
                {
                        is >> stmp >> this->com_remove;
                }
		else if (s.find("cutoff", 0) != string::npos)
                {
                        is >> stmp >> this->cutoff;
                }
		else if (s.find("boundaryType", 0) != string::npos)
		{
			is >> stmp >> this->boundaryType;
		}
		else if (s.find("box_size_x", 0) != string::npos)
		{
			is >> stmp >> this->box_size_x;
		}
		else if (s.find("box_size_y", 0) != string::npos)
		{
			is >> stmp >> this->box_size_y;
		}
		else if (s.find("box_size_z", 0) != string::npos)
		{
			is >> stmp >> this->box_size_z;
		}
		else if (s.find("ewald_kmax", 0) != string::npos)
		{
			is >> stmp >> this->ewald_kmax;
		}
		else if (s.find("ewald_tolerance", 0) != string::npos)
		{
			is >> stmp >> this->ewald_tolerance;
		}
		else if (s.find("usePME", 0) != string::npos)
		{
			is >> stmp >> this->usePME;
		}
		else if (s.find("pme_grid_x", 0) != string::npos)
		{
			is >> stmp >> this->pme_grid_x;
		}
		else if (s.find("pme_grid_y", 0) != string::npos)
		{
			is >> stmp >> this->pme_grid_y;
		}
		else if (s.find("pme_grid_z", 0) != string::npos)
		{
			is >> stmp >> this->pme_grid_z;
		}
		else if (s.find("pme_spline_order", 0) != string::npos)
		{
			is >> stmp >> this->pme_spline_order;
		}
		else if (s.find("wrapAll", 0) != string::npos)
		{
			is >> stmp >> this->wrapAll;
		}
	}

	// consistency check
	//
	if (this->restartfilename.size() && this->initialTemp > 0)
	{
		cerr << "error: restartfilename and temperature(initialTemp) are mutually exclusive parameters.\n";
		return false;
	}
	if (this->binvelocities.size() && this->initialTemp > 0)
	{
		cerr << "error: binvelocities and temperature(initialTemp) are mutually exclusive parameters.\n";
		return false;
	}

	// set seed value if iseed .eq. -1
	//
	if (this->iseed == -1)
	{
		this->iseed = time(NULL);
	}

	if (this->boundaryType == "PBC" && !this->restartfilename.size())
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
	
	if (this->restartfilename.size())
	{
		cout
		<< "REMARK RESTART CONFIGURATION  " << this->restartfilename
		<< '\n'
		<< "REMARK COORDINATE PDB WILL BE IGNORED."
		<< '\n';
	}

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
		if (this->restartfilename.size())
		{
			cout << "REMARK Periodic boundary is assumed.\n";
			cout << "REMARK Box info is loaded from rst7.\n";
		}
		else {
		cout << "REMARK Periodic boundary condition:\n"
			"REMARK     box_size_x: " << this->box_size_x << '\n' <<
			"REMARK     box_size_y: " << this->box_size_y << '\n' <<
			"REMARK     box_size_z: " << this->box_size_z << '\n';
		}
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
