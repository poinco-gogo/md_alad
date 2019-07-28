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

	ComputeBond              tmp_bnd(opt, psf.bondVector);
	ComputeAngle             tmp_ang(psf.angleVector);
	ComputeLJ                tmp_lj(opt, sys, atomVector);
	ComputeES                tmp_es(opt, sys, atomVector, psf);

	this->vbnd               = tmp_bnd;
	this->vang               = tmp_ang;
	this->vlj                = tmp_lj;
	this->ves                = tmp_es;
}

void Energy::calc_temperature()
{
	temperature = kinetic * 2. / ptr_sys->nfree * INVBOLTZMAN;
}

void Energy::calc_kinetic_energy()
{
	this->kinetic = 0;

	for (auto& atom: *ptr_atomVector)
	{
		this->kinetic += atom.mass * atom.velocity.squaredNorm();
	}

	this->kinetic *= 0.5;
}

void Energy::calc_force()
{
	ebond                    = vbnd.compute_force();
	eangle                   = vang.compute_force();
	elj                      = vlj.compute_force();
	ees                      = ves.compute_force();
}

void Energy::zero_force()
{
	for (auto& atom: *ptr_atomVector)
		atom.force = V3ZERO;
}
