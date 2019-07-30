#ifndef ___CLASS_COMPUTE_ES
#define ___CLASS_COMPUTE_ES

#define FC_SYMBOL 1
#include "helpme_standalone.h"
#include <Eigen/Core>
#include "Atom.hpp"
#include "Option.hpp"
#include "System.hpp"
#include "PSF.hpp"
#include "Lattice.hpp"
#include "common.hpp"

class ComputeES
{
	private:
	std::vector<Eigen::Vector3d> g;
	PSF* ptr_psf;
	std::vector<Atom>* ptr_atomVector;
	std::string boundaryType;
	Lattice* ptr_lattice;
	bool usePME;
	int pme_grid_x, pme_grid_y, pme_grid_z, pme_spline_order;
	double cutoff, ewcoef, ewald_kmax, cutoff2, ewcoef2, ewald_tolerance;
	double ewald_self_energy;

	public:
	double tensor[6];
	ComputeES(){}
	ComputeES(const Option& opt, System& sys, std::vector<Atom>& atomVector, PSF& psf);
	double compute_force();

	private:

	void make_reciprocal_vectors();
	void show_simulation_info_ewald();
	double calc_ewcoef();
	double compute_ewald_force();
	double compute_nobc_force();
	double compute_pair_energy(Atom& at1, Atom& at2);
	Eigen::Vector3d compute_pair_force(Atom& at1, Atom& at2);

	double calc_ewald_self();
	double calc_ewald_recip_direct();
	double calc_ewald_recip_pme();
};
#endif
