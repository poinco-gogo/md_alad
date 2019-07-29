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
        int nstep, firsttimestep, DCDFreq;

	std::string rigidBonds, rigidIndexes;
        double rigidTolerance;
        int    rigidIterations;

	std::string boundaryType;
	bool wrapAll;
	double box_size_x, box_size_y, box_size_z;

	int iseed;

	Lattice lattice;

	std::vector<int> soluteIndex, waterIndex;

	System(const Option& opt);

	private:
	void show_simulation_info();
	void make_lattice(double x, double y, double z);
	void make_solute_index(std::vector<Atom>& atomVector);
	void make_water_index(std::vector<Atom>& atomVector);

	public:

	int nfree;
};
#endif
