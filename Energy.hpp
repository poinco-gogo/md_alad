#ifndef ___CLASS_ENERGY
#define ___CLASS_ENERGY

#include <vector>
#include <string>
#include "Atom.hpp"
#include "Lattice.hpp"
#include "Option.hpp"
#include "System.hpp"

class Energy
{
	public:

        int outputEnergies;

	double cutoff;

	std::string usePME;
	int ewald_kmax, pme_grid_x, pme_grid_y, pme_grid_z, pme_spline_order;
	double ewald_tolerance;
	const double ewcoeff = 0.39467;

	Energy(const Option& opt, std::vector<Atom>* ptr_atomVector, System* ptr_sys);

	double es, lj, kinetic;

	void calc_kinetic_energy();
	void calc_potential_energy();
	void calc_force();

	private:
	void show_simulation_info();
	void make_reciprocal_vectors();
	void make_lj_pair_list();
	void make_el_pair_list();

	System* ptr_sys;
	std::vector<Atom>* ptr_atomVector;
	std::vector<int> lj_pair_list, el_pair_list;
	std::vector<Eigen::Vector3d> g;
};
#endif
