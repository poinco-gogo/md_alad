#include <iostream>
#include <iomanip>
#include <cmath>
#include "common.hpp"
#include "ComputeBond.hpp"
using namespace std;

ComputeBond::ComputeBond(const Option& opt, vector<Bond>& bondVector)
{
	this -> ptr_bondVector  = &bondVector;

	this->dt                = opt.dt_fs * 0.001 * PS2ASU;

	this -> rigidBonds      = opt.rigidBonds;
	this -> rigidTolerance  = opt.rigidTolerance;
	this -> rigidIterations = opt.rigidIterations;

	cout << "REMARK RIGID BONDS : ";
	if (rigidBonds == "yes")
	{
		cout << "ALL\n";
		cout << "REMARK         ERROR TOLERANCE : " << rigidTolerance << '\n';
		cout << "REMARK          MAX ITERATIONS : " << rigidIterations << '\n';
		make_shake_bond_index();
	}
	else
	{
		cout << "NONE.\n";
	}
}

double ComputeBond::compute_force()
{
	double sum_energy = 0;

	for (Bond& bond: *ptr_bondVector)
	{
		// skip SHAKE pairs
		//if (rigidBonds == "yes" && bond.is_XH) continue;

		bond.calc_force();
		bond.ptr_atom1->force += bond.force;
		bond.ptr_atom2->force -= bond.force;
		sum_energy += bond.energy;
	}
	return sum_energy;
}

double ComputeBond::do_shake_step()
{
	// loop over all SHAKE bonds.
	for (int i = 0; i < shake_bond_index.size(); i++)
	{
		Bond& bond = ptr_bondVector->at( shake_bond_index[i] );

		Atom* p1 = bond.ptr_atom1;
		Atom* p2 = bond.ptr_atom2;

		Eigen::Vector3d r12new = p1->rnew     - p2->rnew;
		Eigen::Vector3d r12now = p1->position - p2->position;

		double im = p1->invmass + p2->invmass;

		double numer = - ( r12new.squaredNorm() - bond.b0 * bond.b0 );

		double denom = 4. * im * r12new.dot(r12now);

		// Lagrange multiplier
		double dlambda = numer / denom;

		// update the positions
		p1->rnew = p1->rnew + p1->invmass * dlambda * 2. *   r12now;
		p2->rnew = p2->rnew + p2->invmass * dlambda * 2. * (-r12now);
	}

	// error check
	double error = -1;
	for (int i = 0; i < shake_bond_index.size(); i++)
	{
		Bond& bond = ptr_bondVector->at( shake_bond_index[i] );

		Atom* p1 = bond.ptr_atom1;
		Atom* p2 = bond.ptr_atom2;

		Eigen::Vector3d r12new = p1->rnew     - p2->rnew;

		error = max(abs(r12new.squaredNorm()  - bond.b0 * bond.b0), error);
	}

	return error;
}

double ComputeBond::do_rattle_step1()
{
	double invdt = 1. / this->dt;
	double sqinvdt = invdt * invdt;

	// loop over all SHAKE bonds.
	for (int i = 0; i < shake_bond_index.size(); i++)
	{
		Bond& bond = ptr_bondVector->at( shake_bond_index[i] );

		Atom* p1 = bond.ptr_atom1;
		Atom* p2 = bond.ptr_atom2;

		Eigen::Vector3d r12new = p1->rnew     - p2->rnew;
		Eigen::Vector3d r12now = p1->position - p2->position;

		double im = p1->invmass + p2->invmass;

		double numer = - ( r12new.squaredNorm() - bond.b0 * bond.b0 );

		double denom = 4. * im * r12new.dot( r12now );

		// Lagrange multiplier
		double dlambda = numer / denom;

		// multiplier for velocities
		double vlambda = dlambda * invdt;

		// update the positions
		p1->rnew = p1->rnew + p1->invmass * dlambda * 2. *   r12now;
		p2->rnew = p2->rnew + p2->invmass * dlambda * 2. * (-r12now);

		p1->vnew = p1->vnew + p1->invmass * vlambda * 2. *   r12now;
		p2->vnew = p2->vnew + p2->invmass * vlambda * 2. * (-r12now);

		// store and accumulate constraint force
		p1->constfrc += 2.0 * sqinvdt * dlambda * 2. *   r12now;
		p2->constfrc += 2.0 * sqinvdt * dlambda * 2. * (-r12now);
	}

	// error check
	double error = -1;
	for (int i = 0; i < shake_bond_index.size(); i++)
	{
		Bond& bond = ptr_bondVector->at( shake_bond_index[i] );

		Atom* p1 = bond.ptr_atom1;
		Atom* p2 = bond.ptr_atom2;

		Eigen::Vector3d r12new = p1->rnew     - p2->rnew;

		error = max(abs(r12new.squaredNorm()  - bond.b0 * bond.b0), error);
	}

	return error;
}

double ComputeBond::do_rattle_step2()
{
	// loop over all SHAKE bonds.
	for (int i = 0; i < shake_bond_index.size(); i++)
	{
		Bond& bond = ptr_bondVector->at( shake_bond_index[i] );

		Atom* p1 = bond.ptr_atom1;
		Atom* p2 = bond.ptr_atom2;

		// position was already corrected in rattle_vv1()
		Eigen::Vector3d r12now = p1->position - p2->position;

		Eigen::Vector3d v12    = p1->vnew     - p2->vnew;

		double im = p1->invmass + p2->invmass;

		double numer = - r12now.dot( v12 );

		double denom = im * r12now.squaredNorm();

		// Lagrange multiplier
		double dlambda = numer / denom;

		// update the velocities
		p1->vnew = p1->vnew + p1->invmass * dlambda * r12now;
		p2->vnew = p2->vnew + p2->invmass * dlambda * (-r12now);
	}

	// error check
	double error = -1;
	for (int i = 0; i < shake_bond_index.size(); i++)
	{
		Bond& bond = ptr_bondVector->at( shake_bond_index[i] );

		Atom* p1 = bond.ptr_atom1;
		Atom* p2 = bond.ptr_atom2;

		Eigen::Vector3d r12now = p1->position     - p2->position;
		Eigen::Vector3d v12    = p1->vnew     - p2->vnew;

		error = max(abs(v12.dot(r12now)), error);
	}

	return error;
}

bool ComputeBond::do_shake_loop()
{
	bool converged = false;

	for (int ishake = 0; ishake < rigidIterations; ishake++)
	{
		double error = do_shake_step();

		if (error < this->rigidTolerance)
		{
			converged = true;
			break;
		}
	}

	return converged;
}

bool ComputeBond::do_rattle_loop1()
{
	bool converged = false;

	for (int ishake = 0; ishake < rigidIterations; ishake++)
	{
		double error = do_rattle_step1();

		if (error < this->rigidTolerance)
		{
			converged = true;
			break;
		}
	}

	return converged;
}

bool ComputeBond::do_rattle_loop2()
{
	bool converged = false;

	for (int ishake = 0; ishake < rigidIterations; ishake++)
	{
		double error = do_rattle_step2();

		if (error < this->rigidTolerance)
		{
			converged = true;
			break;
		}
	}

	return converged;
}

void ComputeBond::show_bond(const Bond& bond)
{
	cout 
	<< setw(8) << bond.ptr_atom1->PSFIndex
	<< setw(4) << bond.ptr_atom1->PSFAtomName
	<< setw(8) << bond.ptr_atom2->PSFIndex
	<< setw(4) << bond.ptr_atom2->PSFAtomName
	<< setw(8) << bond.Kb
	<< '\n';
}

void ComputeBond::make_shake_bond_index()
{
	for (int i = 0; i < ptr_bondVector->size(); i++)
	{
		if (ptr_bondVector->at(i).is_XH)
			shake_bond_index.push_back(i);
	}
	cout << "REMARK " << shake_bond_index.size() << " shake bond(s) found.\n";
}
