#include <iostream>
#include "Random.hpp"
using namespace std;

Random::Random(unsigned int iseed)
{
	this->iseed       = iseed;
	this->engine.seed(iseed);

	normal_distribution<double> dtmp(0., 1.);
	this->dist = dtmp;
}

double Random::gaussian()
{
	return dist(engine);
}

Eigen::Vector3d Random::gaussian_vector()
{
	return Eigen::Vector3d( gaussian(), gaussian(), gaussian() );
}
