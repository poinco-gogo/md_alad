#ifndef ___CLASS_COMPUTE_LJ
#define ___CLASS_COMPUTE_LJ

#include <Eigen/Core>
#include "Option.hpp"
#include "System.hpp"
#include "Atom.hpp"
#include "Lattice.hpp"

class ComputeLJ
{
	private:
	std::vector<Atom>* ptr_atomVector;
	std::string boundaryType;
	double cutoff, cutoff2;
	Lattice* ptr_lattice;

	public:
	ComputeLJ() {}
	ComputeLJ(const Option& opt, System& sys, std::vector<Atom>& atomVector);
	double compute_force();
	Eigen::Vector3d compute_pair_force(Atom& at1, Atom& at2);

	private:
};
#endif
