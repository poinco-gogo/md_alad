#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include "Integrator.hpp"
#include "common.hpp"
using namespace std;

Integrator::Integrator(const Option& opt)
{
	this->initialTemp        = opt.initialTemp;

	this->dt_fs              = opt.dt_fs;

	this->langevinTemp       = opt.langevinTemp;
	this->langevinDamping_ps = opt.langevinDamping_ps;
	this->langevin           = opt.langevin;
}
