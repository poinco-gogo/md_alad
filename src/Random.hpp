#ifndef ___CLASS_RANDOM
#define ___CLASS_RANDOM

#include <random>
#include "Eigen/Core"

class Random
{
	private:
	
	unsigned int iseed;
	std::mt19937 engine;
	std::normal_distribution<double> dist;

	public:

	Random(unsigned int iseed);

	double gaussian();
	Eigen::Vector3d gaussian_vector();
};
#endif
