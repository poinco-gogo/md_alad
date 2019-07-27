#include <iostream>
#include <iomanip>
#include "Output.hpp"
using namespace std;

Output::Output(ofstream* ptr_fo, System* ptr_sys, Energy* ptr_ene, vector<Atom>* ptr_atomVector)
{
	this->ptr_atomVector     = ptr_atomVector;
	this->ptr_sys            = ptr_sys;
	this->ptr_ene            = ptr_ene;
	this->ptr_fo             = ptr_fo;
}

void Output::print_energy(int istep)
{
	Energy& ene = *ptr_ene;

	double totpot = ene.ebond + ene.eangle + ene.ees + ene.elj;
	double totene = totpot + ene.kinetic;
	double temp   = ene.kinetic * 2. / ptr_sys->nfree * INVBOLTZMAN;

	cout 
		<< setprecision(4) << fixed
		<< setw(12) << istep
		<< setw(16) << ene.ebond
		<< setw(16) << ene.eangle
		<< setw(16) << 0.
		<< setw(16) << 0.
		<< setw(16) << 0.
		<< setw(16) << ene.elj
		<< setw(16) << ene.ees
		<< setw(16) << totpot
		<< setw(16) << ene.kinetic
		<< setw(16) << totene
		<< setw(16) << temp
		<< '\n';
}

void Output::output_xyz()
{
	ofstream& fo = *ptr_fo;
	static const double boxsize = ptr_sys->box_size_x;
	fo << ptr_atomVector->size() << '\n' << '\n';
	fo << setprecision(6) << fixed;
	for (int i = 0; i < ptr_atomVector->size() / 3; i++) 
	{
		int j = 3 * i;
		Atom& at1 = ptr_atomVector->at(j);
		Atom& at2 = ptr_atomVector->at(j + 1);
		Atom& at3 = ptr_atomVector->at(j + 2);
		Eigen::Vector3d p1 = at1.position;
		Eigen::Vector3d p2 = at2.position;
		Eigen::Vector3d p3 = at3.position;
		Eigen::Vector3d vcom = (p1 + p2 + p3) / 3.;

		double d = floor(0.5 + vcom.x() / boxsize);
		p1.x() -= boxsize * d;
		p2.x() -= boxsize * d;
		p3.x() -= boxsize * d;
		
		d = floor(0.5 + vcom.y() / boxsize);
		p1.y() -= boxsize * d;
		p2.y() -= boxsize * d;
		p3.y() -= boxsize * d;
		
		d = floor(0.5 + vcom.z() / boxsize);
		p1.z() -= boxsize * d;
		p2.z() -= boxsize * d;
		p3.z() -= boxsize * d;

		fo << " " << at1.PDBAtomName[0] << setw(20) << p1.x() << setw(20) << p1.y() << setw(20) << p1.z() << '\n';
		fo << " " << at2.PDBAtomName[0] << setw(20) << p2.x() << setw(20) << p2.y() << setw(20) << p2.z() << '\n';
		fo << " " << at3.PDBAtomName[0] << setw(20) << p3.x() << setw(20) << p3.y() << setw(20) << p3.z() << '\n';
	}
}
