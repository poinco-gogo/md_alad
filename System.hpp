#ifndef ___CLASS_SYSTEM
#define ___CLASS_SYSTEM

#include <vector>
#include <string>
#include "Atom.h"
#include "Lattice.hpp"

class System
{
	public:
	std::string filename, restartfilename;
	std::string structure, coordinates, parameters, topparfile,
		outputname, com_remove, boundaryType,
		bincoordinates, binvelocities, wrapAll,
		usePME;
	std::string rigidBonds, rigidIndexes;
        int nstep, firsttimestep,
	    DCDFreq, outputEnergies;
        double initialTemp, cutoff,
	       box_size_x,
	       box_size_y,
	       box_size_z;
        double rigidTolerance;
        int    rigidIterations;
	Lattice lattice;
	int ewald_kmax, pme_grid_x, pme_grid_y, pme_grid_z, pme_spline_order;
	double ewald_tolerance;
	std::vector<int> soluteIndex, waterIndex;

	int iseed;

	System(std::string filename);
	bool load_config();
	void show_simulation_info();
	void make_lattice(double x, double y, double z);
	void make_solute_index(std::vector<Atom>& atomVector);
	void make_water_index(std::vector<Atom>& atomVector);
	
	private:
	void reset();
};
#endif
