#ifndef ___CLASS_ENERGY
#define ___CLASS_ENERGY

#include <vector>
#include <string>
#include "Atom.hpp"
#include "Lattice.hpp"
#include "Option.hpp"

class Energy
{
	public:

        int outputEnergies;

	double cutoff;

	std::string usePME;
	int ewald_kmax, pme_grid_x, pme_grid_y, pme_grid_z, pme_spline_order;
	double ewald_tolerance;
	const double ewcoeff = 0.39467;

	Energy(const Option& opt);

	double es, lj, kinetic;

	void calc_kinetic(std::vector<Atom>& atomVector);

	private:
	void show_simulation_info();
};
#endif
