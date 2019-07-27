#ifndef ___CLASS_COMPUTE_ANGLE
#define ___CLASS_COMPUTE_ANGLE

#include <vector>
#include "Angle.hpp"

class ComputeAngle
{
	private:

	std::vector<Angle>* ptr_angleVector;

	double sum_vangle;
	double sum_vub;
	
	public:

	// constructor
	ComputeAngle(){};
	ComputeAngle(std::vector<Angle>& angleVector);
	void reset();
	double compute_force();

};
#endif
