#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include "Integrator.hpp"
#include "common.hpp"
using namespace std;

Integrator::Integrator(const Option& opt, Energy* ptr_ene, Output* ptr_out, vector<Atom>* ptr_atomVector, mt19937* ptr_engine)
{
	this->integrator         = opt.integrator;

	this->print_energy_step  = opt.outputEnergies;
	this->print_trj_step     = opt.DCDFreq;

	this->initialTemp        = opt.initialTemp;

	this->dt_fs              = opt.dt_fs;

	this->rigidBonds         = opt.rigidBonds;

	this->langevinTemp       = opt.langevinTemp;
	this->langevinDamping_ps = opt.langevinDamping_ps;
	this->langevin           = opt.langevin;

	this->ptr_atomVector     = ptr_atomVector;

	this->ptr_ene            = ptr_ene;
	this->ptr_out            = ptr_out;

	this->ptr_engine         = ptr_engine;

	set_derived_values();

	if (this->langevin) set_langevin_parameters();
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
	normal_distribution<double> dist(0., 1.);

	const double kbT = BOLTZMAN * langevinTemp;

	for (auto& at: *ptr_atomVector)
	{
		at.R = sqrt( 2. * kbT * gamma * at.mass / dt );

		// initial random forces for velocity velret
		at.random_f.x() = at.R * dist(*ptr_engine);
		at.random_f.y() = at.R * dist(*ptr_engine);
		at.random_f.z() = at.R * dist(*ptr_engine);
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

	ptr_ene->calc_kinetic_energy();
	ptr_ene->calc_temperature();

	double factor = sqrt( langevinTemp / ptr_ene->temperature );

	scale_velocity( factor );
}

void Integrator::scale_velocity(const double factor)
{
	for (auto& at: *ptr_atomVector)
		at.velocity *= factor;
}

void Integrator::do_md_loop(const int nstep)
{
	if (this->integrator == "POSI")
		run_position_velret(nstep);
	else if (this->integrator == "VVER")
		run_velocity_velret(nstep);
}

void Integrator::run_position_velret(const int nstep)
{
	Energy& ene = *ptr_ene;
	Output& out = *ptr_out;

	initial_posi_velret();

	for (int istep = 1; istep <= nstep; istep++)
	{
		ene.zero_force();
		ene.calc_force();

		position_velret_integrate();

		if (rigidBonds)
			if (!ene.vbnd.do_shake_loop())
				die("error: shake does not converged!");

		position_velret_velocity();

		ene.calc_kinetic_energy();

		if (istep % print_energy_step== 0)
			out.print_energy(istep);

		if (istep % print_trj_step== 0)
			out.output_dcd();

		for (auto& at: *ptr_atomVector)
		{
			at.rold = at.position;
			at.position = at.rnew;
		}
	}
}

void Integrator::initial_posi_velret()
{
	for (auto& at: *ptr_atomVector)
	{
		at.rnew = at.position + dt * at.velocity
			+ dtdt_div2 * at.invmass * at.force;
	}

	if (rigidBonds)
		if (!ptr_ene->vbnd.do_shake_loop())
			die("error: shake does not converged!");

	for (auto& at: *ptr_atomVector)
	{
		at.rold = at.position;
		at.position = at.rnew;
	}
}

void Integrator::position_velret_integrate()
{
	mt19937& engine = *ptr_engine;
	normal_distribution<double> dist(0., 1.);

	for (auto& at: *ptr_atomVector)
	{
		Eigen::Vector3d
		noise( dist(engine), dist(engine), dist(engine));
		at.rnew = 2. * at.position - A * at.rold + dt * dt * at.invmass * (at.force + at.R * noise);
		at.rnew = at.rnew * inB;
	}
}

void Integrator::position_velret_velocity()
{
	for (auto& at: *ptr_atomVector)
		at.velocity = 0.5 / dt * (at.rnew - at.rold);
}

void Integrator::run_velocity_velret(const int nstep)
{
	Energy& ene = *ptr_ene;
	Output& out = *ptr_out;

	for (int istep = 1; istep <= nstep; istep++)
	{
		velocity_velret_integrate_step1();

		ene.zero_force();
		ene.calc_force();

		velocity_velret_integrate_step2();

		ene.calc_kinetic_energy();

		if (istep % print_energy_step== 0)
			out.print_energy(istep);

		if (istep % print_trj_step== 0)
			out.output_dcd();
	}
}

void Integrator::velocity_velret_integrate_step1()
{
	for (auto& at: *ptr_atomVector)
	{
		at.vnew = A * at.velocity
			+ at.invmass * dt * 0.5 * (at.force + at.random_f);

		at.rnew = at.position + dt * at.vnew;
	}

	if (rigidBonds && !ptr_ene->vbnd.do_rattle_loop1())
		die("error: rattle does not converged!");

	for (auto& at: *ptr_atomVector)
		at.position = at.rnew;
}

void Integrator::velocity_velret_integrate_step2()
{
	normal_distribution<double> dist(0., 1.);

	for (auto& at: *ptr_atomVector)
	{
		at.random_f.x() = at.R * dist(*ptr_engine);
		at.random_f.y() = at.R * dist(*ptr_engine);
		at.random_f.z() = at.R * dist(*ptr_engine);

		at.vnew += at.invmass * dt * 0.5 * (at.force + at.random_f);

		at.vnew *= inB;
	}

	if (rigidBonds && !ptr_ene->vbnd.do_rattle_loop2())
		die("error: rattle does not converged!");

	for (auto& at: *ptr_atomVector)
		at.velocity = at.vnew;
}
