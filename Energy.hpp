#ifndef ___CLASS_ENERGY
#define ___CLASS_ENERGY

#include <vector>
#include <string>
#include "Atom.hpp"
#include "PSF.hpp"
#include "Lattice.hpp"
#include "Option.hpp"
#include "System.hpp"
#include "ComputeBond.hpp"
#include "ComputeAngle.hpp"
#include "ComputeLJ.hpp"
#include "ComputeES.hpp"

class Energy
{
	public:

        int outputEnergies;

	Energy(const Option& opt, System& sys, std::vector<Atom>& atomVector, PSF& psf);

	double temperature;
	double ebond, eangle, ees, elj, kinetic;

	ComputeBond              vbnd;
	ComputeAngle             vang;
	ComputeLJ                vlj;
	ComputeES                ves;

	void calc_temperature();
	void calc_kinetic_energy();
	void calc_force();
	void zero_force();

	private:

	System* ptr_sys;
	std::vector<Atom>* ptr_atomVector;
};
#endif
