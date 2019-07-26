#ifndef ___CLASS_INTEGRATOR
#define ___CLASS_INTEGRATOR

#include <vector>
#include <string>
#include <random>
#include "Atom.hpp"
#include "Option.hpp"

class Integrator
{
	public:

        double initialTemp;

	double langevinTemp, langevinDamping_ps;
	std::string langevin;

	Integrator(const Option& opt, std::vector<Atom>* ptr_atomVector);

	public:

	void set_derived_values();
	void set_langevin_parameters();
	void set_ptr_engine(std::mt19937* ptr_engine);
	void reassign_velocities();
	void initial_posi_velret();

	std::mt19937* ptr_engine;

	std::vector<Atom>* ptr_atomVector;

	double dt_fs, dt_ps, dt, dt_div2, dtdt_div2;
	double gamma_ps, gamma, A, B, inB;
};
#endif
