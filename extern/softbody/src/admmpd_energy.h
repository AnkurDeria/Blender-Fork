// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_ENERGY_H_
#define ADMMPD_ENERGY_H_ 1

#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <vector>

namespace admmpd {

class Lame {
public:
	int m_model; // 0=ARAP
	double m_mu;
	double m_lambda;
	double m_bulk_mod;
	void set_from_youngs_poisson(double youngs, double poisson);
	Lame();
};

class EnergyTerm {
public:

	void signed_svd(
		const Eigen::Matrix<double,3,3> &A, 
		Eigen::Matrix<double,3,3> &U, 
		Eigen::Matrix<double,3,1> &S, 
		Eigen::Matrix<double,3,3> &V);

	// Updates the z and u variables for an element energy.
	void update(
		int index,
		const Lame &lame,
		double rest_volume,
		double weight,
		const Eigen::MatrixXd *x,
		const Eigen::MatrixXd *Dx,
		Eigen::MatrixXd *z,
		Eigen::MatrixXd *u);

	// Initializes tet energy, returns num rows for D
	int init_tet(
		int index,
		const Lame &lame,
		const Eigen::RowVector4i &prim,
		const Eigen::MatrixXd *x,
		double &volume,
		double &weight,
		std::vector< Eigen::Triplet<double> > &triplets);

};

} // end namespace admmpd

#endif // ADMMPD_ENERGY_H_




