#ifndef ___CLASS_INTEGRATOR
#define ___CLASS_INTEGRATOR

#include <vector>
#include <string>
#include <random>
#include "Atom.hpp"
#include "Option.hpp"
#include "Output.hpp"
#include "Energy.hpp"

class Integrator
{
	private:

	std::mt19937* ptr_engine;

	std::vector<Atom>* ptr_atomVector;

	Energy* ptr_ene;

	Output* ptr_out;

	std::string integrator;

	int print_energy_step, print_trj_step;

        double initialTemp;

	double dt_fs, dt_ps, dt, dt_div2, dtdt_div2;
	double gamma_ps, gamma, A, B, inB;

	std::string rigidBonds;

	double langevinTemp, langevinDamping_ps;
	std::string langevin;

	public:

	Integrator(const Option& opt, Energy* ptr_ene, Output* ptr_out, std::vector<Atom>* ptr_atomVector);
	void set_ptr_engine(std::mt19937* ptr_engine);
	void reassign_velocities();
	void do_md_loop(const int nstep);

	private:

	void set_derived_values();
	void set_langevin_parameters();

	void run_position_velret(const int nstep);
	void initial_posi_velret();
	void position_velret_integrate();
	void position_velret_velocity();

	void run_velocity_velret(const int nstep);
};
#endif
