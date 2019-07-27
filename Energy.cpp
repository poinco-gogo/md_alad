#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include "Energy.hpp"
#include "common.hpp"
using namespace std;

Energy::Energy(const Option& opt, System& sys, vector<Atom>& atomVector, PSF& psf)
{
	this->ptr_sys            = &sys;
	this->ptr_atomVector     = &atomVector;

	this->outputEnergies     = opt.outputEnergies;

	ComputeLJ tmp_lj(opt, sys, atomVector);
	ComputeES tmp_es(opt, sys, atomVector, psf);

	this->vlj = tmp_lj;
	this->ves = tmp_es;
}

void Energy::calc_kinetic_energy()
{
	this->kinetic = 0;

	for (auto& atom: *ptr_atomVector)
	{
		this->kinetic += atom.mass * atom.vnew.squaredNorm();
	}

	this->kinetic *= 0.5;
}

void Energy::calc_force()
{
	lj = vlj.compute_force();
	es = ves.compute_force();
}

void Energy::zero_force()
{
	for (auto& atom: *ptr_atomVector)
		atom.force = V3ZERO;
}
