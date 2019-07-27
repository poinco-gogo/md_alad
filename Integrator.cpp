#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include "Integrator.hpp"
#include "common.hpp"
using namespace std;

Integrator::Integrator(const Option& opt, Energy* ptr_ene, Output* ptr_out, vector<Atom>* ptr_atomVector)
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

	set_derived_values();

	if (this->langevin == "yes") set_langevin_parameters();

	make_shake_pairs();
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

void Integrator::do_md_loop(const int nstep)
{
	if (this->integrator == "POSI")
		run_position_velret(nstep);
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

		if (rigidBonds == "yes")
			if (!ene.vbnd.do_shake_loop())
				die("error: shake does not converged!");
/*
		for (int ishake = 0; ishake < 100; ishake++)
		{
			if (shake())
				break;
			if (ishake == 99)
			{
				die("error: shake does not converged.");
			}
		}
*/
		position_velret_velocity();

		ene.calc_kinetic_energy();

		if (istep % print_energy_step== 0)
			out.print_energy(istep);

		if (istep % print_trj_step== 0)
			out.output_xyz();

		for (auto& at: *ptr_atomVector)
		{
			at.rold = at.position;
			at.position = at.rnew;
			at.velocity = at.vnew;
		}
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

void Integrator::position_velret_integrate()
{
	mt19937& engine = *ptr_engine;
	normal_distribution<double> dist(0., 1.);

	for (auto& at: *ptr_atomVector)
	{
		Eigen::Vector3d
		noise( dist(engine), dist(engine), dist(engine));
		at.rnew = 2. * at.position - at.rold + gamma * dt_div2 * at.rold + dt * dt * at.invmass * (at.force + at.R * noise);
		at.rnew = at.rnew * inB;
	}
}

void Integrator::position_velret_velocity()
{
	for (auto& at: *ptr_atomVector)
		at.velocity = 0.5 / dt * (at.rnew - at.rold);
}

bool Integrator::shake()
{
	static const double eps = 1e-6;
	static const double eps2 = eps * eps;
	static const double rOH = 0.9572;
	static const double aHOH = 104.52 * DEG2RAD;
	static const double dOH = rOH * rOH;
	static const double rHH = rOH * sin(aHOH/2.) * 2.;
	static const double dHH = rHH * rHH;

	vector<Atom>& atomVector = *ptr_atomVector;

	for (int i = 0; i < shake_list.size() / 2; i++)
	{
		Atom& at1 = atomVector[shake_list[2 * i]];
		Atom& at2 = atomVector[shake_list[2 * i + 1]];
		double gamma;
		if (at1.PDBAtomName[0] == 'O' && at2.PDBAtomName[0] == 'H')
		{
			gamma = (dOH - (at1.rnew - at2.rnew).squaredNorm()) /
				(2.*(1./at1.mass+1./at2.mass)*((at1.position - at2.position).dot(at1.rnew - at2.rnew)));
		}
		else
		{
			gamma = (dHH - (at1.rnew - at2.rnew).squaredNorm()) /
				(2.*(1./at1.mass+1./at2.mass)*((at1.position - at2.position).dot(at1.rnew - at2.rnew)));

		}
		at1.rnew = at1.rnew + gamma * (at1.position - at2.position) / at1.mass;
		at2.rnew = at2.rnew + gamma * (at2.position - at1.position) / at2.mass;
	}

	for (int i = 0; i < shake_list.size() / 2; i++)
	{
		Atom& at1 = atomVector[shake_list[2 * i]];
		Atom& at2 = atomVector[shake_list[2 * i + 1]];
		double r = (at1.rnew - at2.rnew).norm();
		double error;
		if (at1.PDBAtomName[0]  == 'O' && at2.PDBAtomName[0] == 'H')
		{
			error = abs(r - rOH);
		}
		else
		{
			error = abs(r - rHH);
		}
		if (error > eps) return false;
	}

	return true;
}

void Integrator::make_shake_pairs()
{
	for (int i = 0; i < ptr_atomVector->size() / 3; i++)
	{
		int j = 3 * i;
		shake_list.push_back(j);
		shake_list.push_back(j+1);
		shake_list.push_back(j);
		shake_list.push_back(j+2);
		shake_list.push_back(j+1);
		shake_list.push_back(j+2);
	}
}
