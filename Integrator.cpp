#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include "Integrator.hpp"
#include "common.hpp"
using namespace std;

Integrator::Integrator(const Option& opt, vector<Atom>* ptr_atomVector)
{
	this->initialTemp        = opt.initialTemp;

	this->dt_fs              = opt.dt_fs;

	this->langevinTemp       = opt.langevinTemp;
	this->langevinDamping_ps = opt.langevinDamping_ps;
	this->langevin           = opt.langevin;

	this->ptr_atomVector     = ptr_atomVector;

	set_derived_values();

	if (this->langevin == "yes") set_langevin_parameters();
}

void Integrator::set_derived_values()
{
	this->dt_ps = dt_fs * 0.001;
	this->dt    = dt_ps * PS2ASU;
	this->dt_div2 = dt / 2.;
	this->dtdt_div2 = dt * dt / 2.;

	this->gamma_ps = this->langevinDamping_ps;
	this->gamma    = this->gamma_ps / PS2ASU;
	this->A        = 1. - this->gamma * dt * 0.5;
	this->B        = 1. + this->gamma * dt * 0.5;
	this->inB      = 1. / B;
}

void Integrator::set_langevin_parameters()
{
	for (auto& at: *ptr_atomVector)
	{
		at.R
		= sqrt( 2. * BOLTZMAN * langevinTemp * gamma * at.mass / dt );
	}
}

void Integrator::set_ptr_engine(mt19937* ptr_engine)
{
	this->ptr_engine = ptr_engine;
}

void Integrator::reassign_velocities()
{
	normal_distribution<double> dist(0., 1);

	double kbT = BOLTZMAN * langevinTemp;

	for (auto& at: *ptr_atomVector)
	{
		double kbT_imass = kbT * at.invmass;

		at.velocity.x() = dist( *ptr_engine ) * sqrt( kbT_imass );
		at.velocity.y() = dist( *ptr_engine ) * sqrt( kbT_imass );
		at.velocity.z() = dist( *ptr_engine ) * sqrt( kbT_imass );
	}
}

void Integrator::initial_posi_velret()
{
	for (auto& at: *ptr_atomVector)
	{
		at.fold = at.force;
		at.rold = at.position;
		at.position = at.rold + dt * at.velocity
			+ dtdt_div2 * at.invmass * at.fold;
	}
}
