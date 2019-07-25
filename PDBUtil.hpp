#ifndef ___CLASS_PDBUTIL
#define ___CLASS_PDBUTIL
#include <fstream>
#include <string>
#include "Atom.h"
class PDBUtil
{
	public:
	
	// constructor
	PDBUtil();

	// member funtion
	void writePDBline(const std::ofstream& fo, const Atom& atom);
	void  readPDBline(const  std::string line, const Atom& atom);
};
#endif
