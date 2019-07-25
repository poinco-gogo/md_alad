#ifndef ___CLASS_PSF
#define ___CLASS_PSF

#include <Eigen/Core>
#include <string>
#include <fstream>
#include <vector>
#include "Atom.hpp"

class PSF 
{
	private:
	int natom;
	std::string filename;
	std::vector<Atom>* ptr_atomVector;
	
	public:
	std::vector<Atom>      atomVector;
	
	// constructor
	PSF(std::string filename);
	PSF(std::string filename, std::vector<Atom>* ptr_atomVector);

	// member function
	public:

	void showAtom(const Atom& atom);
	void showResidue(int resID);

	void writePDB(std::string filename, std::string header);
	void writePDB(std::string filename, std::string header, const std::vector<int>& indexVector);
	void writePDB(std::string filename, std::vector<Atom>& atomVector, std::string header);

	// for internal use
	private:
	void read_PSF();
	void get_atom_list(const int natom, std::ifstream& fi);
};
#endif
