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
	std::ofstream* ptr_fo;

	public:
	
	Output(std::ofstream* ptr_fo, System* ptr_sys, Energy* ptr_ene, std::vector<Atom>* ptr_atomVector);
	
	void print_energy(int nstep);
	void output_xyz();

};
#endif
