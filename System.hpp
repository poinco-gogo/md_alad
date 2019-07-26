#ifndef ___CLASS_SYSTEM
#define ___CLASS_SYSTEM

#include <vector>
#include <string>
#include "Atom.hpp"
#include "Lattice.hpp"
#include "Option.hpp"

class System
{
	public:

	std::string structure, coordinates, parameters, topparfile,
		outputname,
		bincoordinates, binvelocities;
        int nstep, firsttimestep, DCDFreq, outputEnergies;

        double initialTemp;

	std::string rigidBonds, rigidIndexes;
        double rigidTolerance;
        int    rigidIterations;

	std::string boundaryType, wrapAll;
	double box_size_x, box_size_y, box_size_z;

	double cutoff;

	std::string usePME;
	int ewald_kmax, pme_grid_x, pme_grid_y, pme_grid_z, pme_spline_order;
	double ewald_tolerance;
	const double ewcoeff = 0.39467;

	int iseed;

	double dt_fs;

	double langevinTemp, langevinDamping_ps;
	std::string langevin;

	Lattice lattice;

	std::vector<int> soluteIndex, waterIndex;

	System(const Option& opt);

	private:
	void show_simulation_info();
	void make_lattice(double x, double y, double z);
	void make_solute_index(std::vector<Atom>& atomVector);
	void make_water_index(std::vector<Atom>& atomVector);
};
#endif
