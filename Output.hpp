#ifndef ___CLASS_OUTPUT
#define ___CLASS_OUTPUT

#include <vector>
#include <fstream>
#include "Option.hpp"
#include "System.hpp"
#include "Energy.hpp"
#include "Atom.hpp"

class Output
{
	private:

	std::vector<Atom>* ptr_atomVector;
	System* ptr_sys;
	Energy* ptr_ene;

	std::ofstream fo;

	std::string outputname;

	public:
	
	Output(const Option& opt, System* ptr_sys, Energy* ptr_ene, std::vector<Atom>* ptr_atomVector);
	
	void print_energy(int nstep);
	void output_xyz();

};
#endif
