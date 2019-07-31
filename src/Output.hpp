#ifndef ___CLASS_OUTPUT
#define ___CLASS_OUTPUT

#include <vector>
#include <fstream>
#include "Option.hpp"
#include "System.hpp"
#include "Energy.hpp"
#include "Atom.hpp"
#include "uiuc/dcdplugin.h"

class Output
{
	private:

	std::vector<Atom>* ptr_atomVector;
	System* ptr_sys;
	Energy* ptr_ene;
	dcdhandle* dcd;
	molfile_timestep_t ts;

	std::ofstream fo;

	std::string outputname;

	public:
	
	Output(const Option& opt, System* ptr_sys, Energy* ptr_ene, std::vector<Atom>* ptr_atomVector);
	
	void print_energy(int nstep);
	void output_xyz();
	void output_dcd();
	bool output_namdbin(std::string type);
	void close_dcd();

	private:

	void open_dcd_write();
};
#endif
