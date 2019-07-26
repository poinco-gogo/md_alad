#ifndef ___CLASS_INTEGRATOR
#define ___CLASS_INTEGRATOR

#include <vector>
#include <string>
#include "Option.hpp"

class Integrator
{
	public:

        double initialTemp;

	double dt_fs;

	double langevinTemp, langevinDamping_ps;
	std::string langevin;

	Integrator(const Option& opt);
};
#endif
