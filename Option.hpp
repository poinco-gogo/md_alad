#ifndef ___CLASS_OPTION
#define ___CLASS_OPTION

#include <string>

class Option
{
	public:
	
	Option(std::string filename);
	bool load_config();
	void reset();

	std::string filename;
	std::string structure, coordinates, parameters, topparfile,
		outputname,
		bincoordinates, binvelocities;
        int nstep, firsttimestep, DCDFreq, outputEnergies;

        double initialTemp;

	std::string integrator;

	std::string rigidBonds, rigidIndexes;
        double rigidTolerance;
        int    rigidIterations;

	std::string boundaryType, wrapAll;
	double box_size_x, box_size_y, box_size_z;

	double cutoff;

	std::string usePME;
	int ewald_kmax, pme_grid_x, pme_grid_y, pme_grid_z, pme_spline_order;
	double ewald_tolerance;

	int iseed;

	double dt_fs;

	double langevinTemp, langevinDamping_ps;
	std::string langevin;
};
#endif
