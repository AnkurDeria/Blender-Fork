


#include "admmpd_energy.h"
#include <iostream>

namespace admmpd {
using namespace Eigen;

Lame::Lame() : m_model(0)
{
	set_from_youngs_poisson(100000,0.299);
}

void Lame::set_from_youngs_poisson(double youngs, double poisson)
{
	m_mu = youngs/(2.0*(1.0+poisson));
	m_lambda = youngs*poisson/((1.0+poisson)*(1.0-2.0*poisson));
	m_bulk_mod = m_lambda + (2.0/3.0)*m_mu;
}

void EnergyTerm::signed_svd(
	const Eigen::Matrix<double,3,3>& A, 
	Eigen::Matrix<double,3,3> &U, 
	Eigen::Matrix<double,3,1> &S, 
	Eigen::Matrix<double,3,3> &V)
{
	JacobiSVD<Matrix3d> svd(A, ComputeFullU | ComputeFullV);
	S = svd.singularValues();
	U = svd.matrixU();
	V = svd.matrixV();
	Matrix3d J = Matrix3d::Identity();
	J(2,2) = -1.0;
	if( U.determinant() < 0.0 )
	{
		U = U * J;
		S[2] = -S[2];
	}
	if( V.determinant() < 0.0 )
	{
		Matrix3d Vt = V.transpose();
		Vt = J * Vt;
		V = Vt.transpose();
		S[2] = -S[2];
	}
}

void EnergyTerm::update(
	int index,
	const Lame &lame,
	double rest_volume,
	double weight,
	const Eigen::MatrixXd *x,
	const Eigen::MatrixXd *Dx,
	Eigen::MatrixXd *z,
	Eigen::MatrixXd *u)
{
	Matrix3d Dix = Dx->block<3,3>(index,0);
	Matrix3d ui = u->block<3,3>(index,0);
	Matrix3d zi = Dix + ui;

	Matrix3d U, V;
	Vector3d S;
	signed_svd(zi, U, S, V);
	S = Vector3d::Ones();

	Matrix3d p = U * S.asDiagonal() * V.transpose();
	double k = lame.m_bulk_mod;
	double kv = k * rest_volume;
	double w2 = weight*weight;
	zi = (kv*p + w2*zi) / (w2 + kv);
	ui += (Dix - zi);

	u->block<3,3>(index,0) = ui;
	z->block<3,3>(index,0) = zi;

} // end EnergyTerm::update

int EnergyTerm::init_tet(
	int index,
	const Lame &lame,
	const Eigen::Matrix<int,1,4> &prim,
	const Eigen::MatrixXd *x,
	double &volume,
	double &weight,
	std::vector< Eigen::Triplet<double> > &triplets )
{
	Matrix<double,3,3> edges;
	edges.col(0) = x->row(prim[1]) - x->row(prim[0]);
	edges.col(1) = x->row(prim[2]) - x->row(prim[0]);
	edges.col(2) = x->row(prim[3]) - x->row(prim[0]);
	Matrix<double,3,3> edges_inv = edges.inverse();
	volume = edges.determinant() / 6.0f;
	if( volume < 0 )
		throw std::runtime_error("**Solver::energy_init: Inverted initial tet");
	double k = lame.m_bulk_mod;
std::cout << "IDX: " << index << " bulk mod: " << k << std::endl;
	weight = std::sqrt(k*volume);
	Matrix<double,4,3> S = Matrix<double,4,3>::Zero();
	S(0,0) = -1; S(0,1) = -1; S(0,2) = -1;
	S(1,0) =  1; S(2,1) =  1; S(3,2) =  1;
	Eigen::Matrix<double,4,3> D = S * edges_inv;
	Eigen::Matrix<double,3,4> Dt = D.transpose();
	int rows[3] = { index+0, index+1, index+2 };
	int cols[4] = { prim[0], prim[1], prim[2], prim[3] };
	for( int r=0; r<3; ++r ){
		for( int c=0; c<4; ++c ){
			triplets.emplace_back(rows[r], cols[c], Dt(r,c));
		}
	}
	return 3;
}

} // end namespace mcl
