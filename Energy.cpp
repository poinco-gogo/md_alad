#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include "Energy.hpp"
#include "common.hpp"
using namespace std;

Energy::Energy(const Option& opt)
{
	this->outputEnergies     = opt.outputEnergies;

	this->cutoff             = opt.cutoff;

	this->usePME             = opt.usePME;
	this->ewald_kmax         = opt.ewald_kmax;
	this->pme_grid_x         = opt.pme_grid_x;
	this->pme_grid_y         = opt.pme_grid_y;
	this->pme_grid_z         = opt.pme_grid_z;
	this->pme_spline_order   = opt.pme_spline_order;
	this->ewald_tolerance    = opt.ewald_tolerance;
}

void Energy::show_simulation_info()
{
	cout 
		<< "REMARK CUTOFF                 " << this->cutoff
		<< '\n';
}
