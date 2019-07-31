#ifndef ___CLASS_COMPUTEBOND
#define ___CLASS_COMPUTEBOND

#include <vector>
#include <string>
#include "Bond.hpp"
#include "Option.hpp"
//#include "System.hpp"

class ComputeBond
{
	private:
	bool rigidBonds;
	double  rigidTolerance;
	int     rigidIterations;

	std::vector<int> shake_bond_index;

	std::vector<Bond>* ptr_bondVector;

	public:
	// constructor
	ComputeBond(){};
	ComputeBond(const Option& opt, std::vector<Bond>& bondVector);

	// member 
	double do_shake_step();
	double do_rattle_step1();
	double do_rattle_step2();
	bool do_shake_loop();
	bool do_rattle_loop1();
	bool do_rattle_loop2();

	void show_bond(const Bond& bond);

	double compute_force();

	int _num_shake_bond(){ return shake_bond_index.size(); }

	private:

	double dt;

	void make_shake_bond_index();
};
#endif
