#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <fstream>
#include <sstream>
#include "PSF.hpp"
#include "PDBUtil.hpp"

using namespace std;

PSF::PSF(string filename) 
{
	ptr_atomVector = &atomVector;
	this -> filename = filename;
	read_PSF();
}

PSF::PSF(string filename, vector<Atom>* ptr_atomVector) 
{
	this -> ptr_atomVector = ptr_atomVector;
	this -> filename = filename;
	read_PSF();
}

void PSF::read_PSF()
{
	ifstream fi(filename.c_str());

	if ( !fi )
	{
		cerr << "エラー：ファイル" << filename << "を開けません。\n";

		return;
	}

	int loc;
	
	while ( getline(fi, s) )
	{
		loc = s.find("NATOM", 0);
		if (loc != string::npos) 
		{
			istringstream is;
			is.str(s.substr(0, 8)); // number of atoms in system
			is >> natom;
			get_atom_list(natom, fi); 
		}
	}

	fi.close();
}

void PSF::get_atom_list(const int natom, const ifstream& fi)
{
	cout << "REMARK TOTAL ATOM : " << natom << '\n';
	ptr_atomVector -> resize(natom);
	string s;
	string stmp;
	int itmp;
	for (int i = 0; i < natom; i++)
	{
		Atom& atom = ptr_atomVector->at(i);
		getline(fi, s);
		istringstream is(s);
		is >> atom.PSFIndex 
		   >> atom.PSFSegmentName
		   >> atom.PSFResID
		   >> atom.PSFResName
		   
		   /*  overwrite PDB value */
		   >> atom.PDBAtomName
		   
		   >> atom.PSFAtomName 
		   >> atom.charge
		   >> atom.mass;

		atom.invmass = 1. / atom.mass;
	}
}

void PSF::writePDB(string filename, string header, const vector<int>& indexVector)
{
	PDBUtil PDBJOBS;
	
	ofstream fo(filename.c_str());

	time_t now = time(NULL);
	struct tm* local = localtime(&now);
	
	fo << "REMARK " << asctime(local);
	fo << "REMARK " << header;
	fo << '\n';
	
	for (int i = 0; i < indexVector.size(); i++)
	{
		int ii = ptr_indexVector -> at(i) - 1;
		PDBJOBS.writePDBline(fo, ptr_atomVector -> at(ii));
	}
}// end of function writePDB()

void PSF::writePDB(string filename, string header)
{
	PDBUtil PDBJOBS;
	
	ofstream fo(filename.c_str());

	time_t now = time(NULL);
	struct tm* local = localtime(&now);
	
	fo << "REMARK " << asctime(local);
	fo << "REMARK " << header;
	fo << '\n';
	
	for (int i = 0; i < ptr_atomVector -> size(); i++)
	{
		PDBJOBS.writePDBline(fo, ptr_atomVector -> at(i));
	}
}// end of function writePDB()

void PSF::writePDB(string filename, const vector<Atom>& atomVector, string header)
{
	PDBUtil PDBJOBS;
	
	ofstream fo(filename.c_str());

	time_t now = time(NULL);
	struct tm* local = localtime(&now);
	
	fo << "REMARK  " << asctime(local);
	fo << "REMARK  " << header;
	fo << '\n';
	
	for (int i = 0; i < atomVector.size(); i++)
	{
		Atom& base = ptr_atomVector->at(i);
		Atom& atom = atomVector[i];

		atom.PSFIndex = base.PSFIndex;
		atom.PDBAtomName = base.PDBAtomName;
		atom.PSFResName = base.PSFResName;
		atom.PSFSegmentName = base.PSFSegmentName;
		atom.PSFResID = base.PSFResID;
//		atom.PDBo = base.PDBo;
//		atom.PDBb = base.PDBb;
		atom.PDBo = 0;
		atom.PDBb = 0;

		PDBJOBS.writePDBline(fo, atomVector[i]);
	}
}// end of function writePDB()

void PSF::showAtom(const Atom& atom)
{
	cout
	<< atom.PSFIndex << ' ' << atom.PDBAtomName
	<< " @" << atom.PSFResName << atom.PSFResID 
	<< "  " << atom.charge << 'e' 
	<< "  " << atom.mass   << "amu\n";
}

void PSF::showResidue(int resID)
{
	bool hit = false;
	for (int i = 0; i < ptr_atomVector -> size(); i++)
	{
		if (ptr_atomVector -> at(i).PSFResID == resID)
		{
			cout << ptr_atomVector -> at(i).PSFResName << '\n';
			hit = true;
			break;
		}
	}
	if (!hit)
		cout << "REMARK warning: no residue hit for resID: "
		<< resID << '\n';
}
