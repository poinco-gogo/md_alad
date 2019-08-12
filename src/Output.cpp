#include <iostream>
#include <iomanip>
#include "Output.hpp"
#include "NAMDBin.hpp"
using namespace std;

Output::Output(const Option& opt, System* ptr_sys, Energy* ptr_ene, vector<Atom>* ptr_atomVector)
{
	this->ptr_atomVector     = ptr_atomVector;
	this->ptr_sys            = ptr_sys;
	this->ptr_ene            = ptr_ene;
	this->outputname         = opt.outputname;

	//this->fo.open(outputname + ".xyz");
	this->open_dcd_write();
}

void Output::open_dcd_write()
{
	char ctmp;
	int with_unitcell = 0;

	if (ptr_sys->boundaryType == "PBC") with_unitcell = 1;

	this->dcd = ::open_dcd_write(
			(outputname + ".dcd").c_str(),
			&ctmp,
			ptr_atomVector->size(),
			with_unitcell );

	ts.coords = (float *)malloc(dcd->natoms * 3 * sizeof(float));
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
		<< setw(16) << ene.edihed
		<< setw(16) << ene.ecmap
		<< setw(16) << ene.eimprop
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

void Output::output_dcd()
{
	for (int i = 0; i < ptr_atomVector->size(); i++)
	{
		Atom& at = ptr_atomVector->at(i);

		ts.coords[3 * i    ] = at.position.x();
		ts.coords[3 * i + 1] = at.position.y();
		ts.coords[3 * i + 2] = at.position.z();
	}

	ts.A = ptr_sys->lattice._x();
	ts.B = ptr_sys->lattice._y();
	ts.C = ptr_sys->lattice._z();

	::write_timestep(dcd, &ts);
}

bool Output::output_namdbin(string type)
{
	NAMDBin job(this->outputname + "." + type, type);

	return job.write_fo(*ptr_atomVector);
}

void Output::close_dcd()
{
	close_file_write( dcd );
	free( ts.coords );
}
