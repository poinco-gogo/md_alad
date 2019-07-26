#ifndef ___CLASS_OUTPUT
#define ___CLASS_OUTPUT

#include <vector>
#include <fstream>
#include "System.hpp"
#include "Energy.hpp"
#include "Atom.hpp"

class Output
{
	private:

	std::vector<Atom>* ptr_atomVector;
	System* ptr_sys;
	Energy* ptr_ene;

	public:
	
	Output(std::vector<Atom>* ptr_atomVector, System* ptr_sys, Energy* ptr_ene);
	
	void print_energy(int nstep);
	void output_xyz(std::ofstream& fo);

};
#endif
