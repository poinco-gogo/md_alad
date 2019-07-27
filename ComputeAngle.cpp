#include <iostream>
#include <iomanip>
#include "ComputeAngle.hpp"
using namespace std;

ComputeAngle::ComputeAngle(vector<Angle>& angleVector)
{
	this -> ptr_angleVector = &angleVector;
	reset();
}

void ComputeAngle::reset()
{
	sum_vangle = 0.0;
	sum_vub    = 0.0;
}

double ComputeAngle::compute_force()
{
	sum_vangle = 0;
	sum_vub    = 0;
	double sum_energy = 0;

	for (Angle& angle: *ptr_angleVector)
	{
		angle.calc_force();

		angle.ptr_atom1->force += angle.force1;
		angle.ptr_atom2->force += angle.force2;
		angle.ptr_atom3->force += angle.force3;

		sum_vangle += angle.energy;
		sum_energy += angle.energy;

		// calculate Urey-Bradley
		if (angle.Kub != 0.0)
		{	
			angle.calc_force_ub();

			angle.ptr_atom1->force += angle.force_ub;
			angle.ptr_atom3->force -= angle.force_ub;
			sum_vub    += angle.energy_ub;
			sum_energy += angle.energy_ub;
		}
	}

	return sum_energy;
}
