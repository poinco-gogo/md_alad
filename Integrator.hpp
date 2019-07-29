#ifndef ___CLASS_INTEGRATOR
#define ___CLASS_INTEGRATOR

#include <vector>
#include <string>
#include <random>
#include "Random.hpp"
#include "Atom.hpp"
#include "Option.hpp"
#include "Output.hpp"
#include "Energy.hpp"

class Integrator
{
	private:

	Random* ptr_random;

	std::vector<Atom>* ptr_atomVector;

	Energy* ptr_ene;

	Output* ptr_out;

	std::string integrator;

	int print_energy_step, print_trj_step;

        double initialTemp;

	double dt_fs, dt_ps, dt, dt_div2, dtdt_div2;
	double gamma_ps, gamma, A, B, inB;

	bool rigidBonds;

	double langevinTemp, langevinDamping_ps;
	bool langevin;

	public:

	Integrator(const Option& opt, Energy* ptr_ene, Output* ptr_out, std::vector<Atom>* ptr_atomVector, Random* ptr_random);

	void reassign_velocities();
	void do_md_loop(const int nstep);

	private:

	void set_derived_values();
	void set_langevin_parameters();

	void scale_velocity(const double factor);

	void run_position_velret(const int nstep);
	void initial_posi_velret();
	void position_velret_integrate();
	void position_velret_velocity();

	void run_velocity_velret(const int nstep);
	void velocity_velret_integrate_step1();
	void velocity_velret_integrate_step2();
};
#endif
