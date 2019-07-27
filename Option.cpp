#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <algorithm>
#include "Option.hpp"
using namespace std;
Option::Option(string filename)
{
	this->filename = filename;
	reset();
}

void Option::reset()
{
	nstep          = 0;
	firsttimestep  = 0;
	DCDFreq        = 1;
	outputEnergies = 1;

	initialTemp    = 0;

	integrator     = "POSI";

	rigidBonds      = false;
	rigidTolerance  = 1e-10;
	rigidIterations = 100;

	boundaryType = "NOBC";
	wrapAll      = "no";
	box_size_x = 0;
	box_size_y = 0;
	box_size_z = 0;

	cutoff = 12.;

	usePME = false;
	ewald_kmax = 10;
	pme_grid_x = 32;
	pme_grid_y = 32;
	pme_grid_z = 32;
	pme_spline_order = 4;
	ewald_tolerance  = 1e-6;

	iseed = -1;

	dt_fs              = 0.5;

	langevinTemp       = 300;
	langevinDamping_ps = 0.5;
	langevin           = false;
}

bool Option::load_config()
{
	ifstream fc(filename.c_str());

	string s;
	while (getline(fc, s))
	{
		istringstream is(s);
		string stmp;
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
                else if (s.find("integrator", 0) != string::npos)
                {
                        is >> stmp >> this->integrator;
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
                        is >> stmp >> stmp;
			this->rigidBonds = stobool(stmp);
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
			is >> stmp >> stmp;
			this->usePME = stobool(stmp);
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
		else if (s.find("timestep", 0) != string::npos)
		{
			is >> stmp >> this->dt_fs;
		}
		else if (s.find("langevinTemp", 0) != string::npos)
		{
			is >> stmp >> this->langevinTemp;
		}
		else if (s.find("langevinDamping", 0) != string::npos)
		{
			is >> stmp >> this->langevinDamping_ps;
		}
		else if (s.find("langevin", 0) != string::npos)
		{
			is >> stmp >> stmp;
			this->langevin = stobool(stmp);
		}
		else
		{
			is >> stmp;
			cerr << "\nunknown config param \n" << stmp
				<< "\nfound.\n\n";
			return false;
		}
	}

	// consistency check
	if (this->binvelocities.size() && this->initialTemp > 0)
	{
		cerr << "\nerror: binvelocities and temperature are mutually exclusive parameters.\n";
		return false;
	}

	// set seed value if iseed .eq. -1
	if (this->iseed == -1)
	{
		this->iseed = time(NULL);
	}
}

bool Option::stobool(string s)
{
	for (auto& c: s) c = toupper(c);

	if (s == "YES" || s == "ON") return true;
	else return false;
}
