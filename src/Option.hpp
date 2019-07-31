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

	bool rigidBonds;
	std::string rigidIndexes;
        double rigidTolerance;
        int    rigidIterations;

	std::string boundaryType;
	bool wrapAll;
	double box_size_x, box_size_y, box_size_z;

	double cutoff;

	bool usePME;
	int ewald_kmax, pme_grid_x, pme_grid_y, pme_grid_z, pme_spline_order;
	double ewald_tolerance;

	unsigned int iseed;

	double dt_fs;

	double langevinTemp, langevinDamping_ps;
	bool langevin;

	double berendsenTemp, berendsenPeriod_ps;
	bool berendsen;

	private:

	bool consistency_check();

	bool stobool(std::string);
};
#endif
