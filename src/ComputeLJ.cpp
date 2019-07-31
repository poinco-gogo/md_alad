#include <iostream>
#include <iomanip>
#include "ComputeLJ.hpp"
using namespace std;

ComputeLJ::ComputeLJ(const Option& opt, System& sys, vector<Atom>& atomVector)
{
	this->ptr_atomVector     = &atomVector;
	this->cutoff             = opt.cutoff;
	this->cutoff2            = cutoff * cutoff;
	this->boundaryType       = sys.boundaryType;
	this->ptr_lattice        = &sys.lattice;
}

double ComputeLJ::compute_force()
{
	Lattice lattice = *ptr_lattice;

	double sum_energy = 0;

	double work[3][3] = {};

	for (int i = 0; i < ptr_atomVector->size(); i++)
	{
		Atom& iat = ptr_atomVector->at(i);

		for (auto& j: iat.ex_pair_list)
		{
			Atom& jat = ptr_atomVector->at(j);

			bool scaled1_4 = false;
			if (jat.checkScaled1_4Pair(iat))
				scaled1_4 = true;
        
			double epsij = 0.;
			double Rmnij = 0.;
			if (scaled1_4)
        		{
                		epsij = sqrt(iat.eps1_4 * jat.eps1_4);
                		Rmnij = iat.Rmin1_4 + jat.Rmin1_4;
        		}
        		else
        		{
                		epsij = sqrt(iat.epsilon * jat.epsilon);
                		Rmnij = iat.Rmin_div2 + jat.Rmin_div2;
        		}

			Eigen::Vector3d r12;

			if (boundaryType == "PBC")
				r12 = lattice.delta(iat.position, jat.position);
			else
				r12 = iat.position - jat.position;

			double rsq = r12.squaredNorm();

			if (rsq > cutoff2) continue;

			double invr2 = 1. / rsq;
			double lj6 = Rmnij * Rmnij;
			lj6 *= invr2;
			lj6 *= lj6 * lj6;
			double lj12 = lj6 * lj6;

			sum_energy += epsij * (lj12 - 2. * lj6);

			Eigen::Vector3d f = epsij * 12. * invr2 * lj6 * (lj6 - 1.) * r12;
			iat.force += f;
			jat.force -= f;
/*
			work[0][0] += f.x() * r12.x();
			work[0][1] += f.x() * r12.y();
			work[0][2] += f.x() * r12.z();
			work[1][0] += f.y() * r12.x();
			work[1][1] += f.y() * r12.y();
			work[1][2] += f.y() * r12.z();
			work[2][0] += f.z() * r12.x();
			work[2][1] += f.z() * r12.y();
			work[2][2] += f.z() * r12.z();
*/
		}
	}

	//cout << setprecision(4) << fixed;
	//cout << "DEBUG: lj: " << work[0][0] << " " << work[1][1] << " " << work[2][2] << '\n';
/*
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			tensor[i][j] += work[i][j];
*/
	return sum_energy;
}

Eigen::Vector3d ComputeLJ::compute_pair_force(Atom& atom1, Atom& atom2)
{
	// return 0 if atom 1 and 2 are 1-2,1-3 pair
        if (atom1.checkExclusionPair(atom2))
        {
                Eigen::Vector3d vtmp(0., 0., 0.);
                return vtmp;
        }
        
	bool scaled1_4 = false;
        // do i have to use scaled 1-4 LJ parameters?
        if (atom1.checkScaled1_4Pair(atom2))
                scaled1_4 = true;

        double epsij = 0.;
	double Rmnij = 0.;
	if (scaled1_4)
        {
                epsij = sqrt(atom1.eps1_4 * atom2.eps1_4);
                Rmnij = atom1.Rmin1_4 + atom2.Rmin1_4;
        }
        else
        {
                epsij = sqrt(atom1.epsilon * atom2.epsilon);
                Rmnij = atom1.Rmin_div2 + atom2.Rmin_div2;
        }
        Eigen::Vector3d r12  = atom1.position - atom2.position;
        double inv2  = 1. / r12.norm();
        inv2 *= inv2;
        double lj6 = Rmnij * Rmnij * inv2;
        lj6 = lj6 * lj6 * lj6;
        double lj12 = lj6 * lj6;
	return  12. * epsij * inv2 * (lj12 - lj6) * r12;
}
