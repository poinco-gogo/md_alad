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

	this->cutoff             = opt.cutoff;

	ComputeBond              tmp_bnd(opt, psf.bondVector);
	ComputeAngle             tmp_ang(psf.angleVector);
	ComputeDihedral          tmp_dih(psf.dihedralVector);
	ComputeImproper          tmp_imp(psf.improperVector);
	ComputeLJ                tmp_lj(opt, sys, atomVector);
	ComputeES                tmp_es(opt, sys, atomVector, psf);

	this->vbnd               = tmp_bnd;
	this->vang               = tmp_ang;
	this->vdih               = tmp_dih;
	this->vimp               = tmp_imp;
	this->vlj                = tmp_lj;
	this->ves                = tmp_es;

	show_simulation_info();

	make_exclusion_pairs();
}

void Energy::show_simulation_info()
{
	cout << "REMARK cutoff[ang.] " << cutoff << '\n';
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
	edihed                   = vdih.compute_force();
	eimprop                  = vimp.compute_force();
	elj                      = vlj.compute_force();
	ees                      = ves.compute_force();
}

void Energy::zero_force()
{
	for (auto& atom: *ptr_atomVector)
		atom.force = V3ZERO;
}

void Energy::make_exclusion_pairs()
{
	for (int i = 0; i < ptr_atomVector->size(); i++)
	{
		Atom& iat = ptr_atomVector->at(i);

		for (int j = i + 1; j < ptr_atomVector->size(); j++)
		{
			Atom& jat = ptr_atomVector->at(j);

			if (jat.checkExclusionPair(iat)) continue;

			iat.ex_pair_list.push_back(j);
		}
	}
}
