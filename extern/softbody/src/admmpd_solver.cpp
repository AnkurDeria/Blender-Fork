


#include "admmpd_solver.h"
#include "admmpd_lattice.h"
#include "admmpd_energy.h"
#include "admmpd_collision.h"

#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <stdio.h>
#include <iostream>

namespace admmpd {
using namespace Eigen;
template <typename T> using RowSparseMatrix = SparseMatrix<T,RowMajor>;

bool Solver::init(
    const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &T,
    const Options *options,
    Data *data)
{
	if (!data || !options)
		throw std::runtime_error("init: data/options null");

	data->x = V;
	data->tets = T;
	compute_matrices(options,data);
	return true;
} // end init

int Solver::solve(
	const Options *options,
	Data *data)
{

	// Init the solve which computes
	// quantaties like M_xbar and makes sure
	// the variables are sized correctly.
	init_solve(options,data);

	int ne = data->rest_volumes.size();
	Lame lame; // TODO lame params as input

	// Begin solver loop
	int iters = 0;
	for (; iters < options->max_admm_iters; ++iters)
	{
		// Local step
		data->Dx.noalias() = data->D * data->x;
		for(int i=0; i<ne; ++i)
		{
			EnergyTerm().update(
				data->indices[i][0],
				lame,
				data->rest_volumes[i],
				data->weights[i],
				&data->x,
				&data->Dx,
				&data->z,
				&data->u );
		}

		// Global step
		update_constraints(options,data);
		data->b.noalias() = data->M_xbar + data->DtW2*(data->z-data->u);
		solve_conjugate_gradients(options,data);

	} // end solver iters

	double dt = options->timestep_s;
	if (dt > 0.0)
		data->v.noalias() = (data->x-data->x_start)*(1.0/dt);

	return iters;
} // end solve

void Solver::init_solve(
	const Options *options,
	Data *data)
{
	int nx = data->x.rows();
	if (data->M_xbar.rows() != nx)
		data->M_xbar.resize(nx,3);

	// velocity and position
	double dt = std::max(0.0, options->timestep_s);
	data->x_start = data->x;
	for (int i=0; i<nx; ++i)
	{
		data->v.row(i) += options->grav;
		data->M_xbar.row(i) =
			data->m[i] * data->x.row(i) +
			dt*data->m[i]*data->v.row(i);
	}

	// ADMM variables
	data->Dx.noalias() = data->D * data->x;
	data->z = data->Dx;
	data->u.setZero();

} // end init solve

void Solver::update_constraints(
	const Options *options,
	Data *data)
{

	std::vector<double> l_coeffs;
	std::vector<Eigen::Triplet<double> > trips_x;
    std::vector<Eigen::Triplet<double> > trips_y;
    std::vector<Eigen::Triplet<double> > trips_z;

	// TODO collision detection
	FloorCollider().jacobian(
		&data->x,
		&trips_x,
		&trips_y,
		&trips_z,
		&l_coeffs);

	// Check number of constraints.
	// If no constraints, clear Jacobian.
	int nx = data->x.rows();
	int nc = l_coeffs.size();
	if (nc==0)
	{
		data->l.setZero();
		for (int i=0; i<3; ++i)
			data->K[i].setZero();

		return;
	}

	// Otherwise update the data.
	data->l = Map<VectorXd>(l_coeffs.data(),nc);
	data->K[0].resize(nc,nx);
	data->K[0].setFromTriplets(trips_x.begin(),trips_x.end());
	data->K[1].resize(nc,nx);
	data->K[1].setFromTriplets(trips_y.begin(),trips_y.end());
	data->K[2].resize(nc,nx);
	data->K[2].setFromTriplets(trips_z.begin(),trips_z.end());

} // end update constraints

void Solver::solve_conjugate_gradients(
	const Options *options,
	Data *data)
{
	// If we don't have any constraints,
	// we don't need to perform CG
	int nnz = std::max(std::max(
		data->K[0].nonZeros(),
		data->K[1].nonZeros()),
		data->K[2].nonZeros());
	if (nnz==0)
	{
		data->x.noalias() = data->ldltA.solve(data->b);
		return;
	}

	// Inner product of matrices interpreted
	// if they were instead vectorized
	auto mat_inner = [](
		const MatrixXd &A,
		const MatrixXd &B )
	{
		double dot = 0.0;
		int nr = std::min(A.rows(), B.rows());
		for( int i=0; i<nr; ++i )
			for(int j=0; j<3; ++j)
				dot += A(i,j)*B(i,j);

		return dot;
	};

	double eps = options->min_res;
	MatrixXd b = data->b;
	int nv = b.rows();
	RowSparseMatrix<double> A[3];
	MatrixXd r(b.rows(),3);
	MatrixXd z(nv,3);
	MatrixXd p(nv,3);
	MatrixXd Ap(nv,3);

	for (int i=0; i<3; ++i)
	{
		RowSparseMatrix<double> Kt = data->K[i].transpose();
		A[i] = data->A + data->spring_k*RowSparseMatrix<double>(Kt*data->K[i]);
		b.col(i) += data->spring_k*Kt*data->l;
		r.col(i) = b.col(i) - A[i]*data->x.col(i);
	}
	z = data->ldltA.solve(r);
	p = z;

	for (int iter=0; iter<options->max_cg_iters; ++iter)
	{
		for( int i=0; i<3; ++i )
			Ap.col(i) = A[i]*p.col(i);

		double p_dot_Ap = mat_inner(p,Ap);
		if( p_dot_Ap==0.0 )
			break;

		double zk_dot_rk = mat_inner(z,r);
		if( zk_dot_rk==0.0 )
			break;

		double alpha = zk_dot_rk / p_dot_Ap;
		data->x += alpha * p;
		r -= alpha * Ap;
		if( r.lpNorm<Infinity>() < eps )
			break;

		z = data->ldltA.solve(r);
		double beta = mat_inner(z,r) / zk_dot_rk;
		p = z + beta*p;
	}
} // end solve conjugate gradients

void Solver::compute_matrices(
	const Options *options,
	Data *data)
{
	// Allocate per-vertex data
	int nx = data->x.rows();
	data->x_start = data->x;
	data->M_xbar.resize(nx,3);
	data->M_xbar.setZero();
	data->Dx.resize(nx,3);
	data->Dx.setZero();
	if (data->v.rows() != nx)
	{
		data->v.resize(nx,3);
		data->v.setZero();
	}
	if (data->m.rows() != nx)
	{ // TODO get from BodyPoint
		data->m.resize(nx);
		data->m.setOnes();
		data->m *= 0.01;
	}

	// Add per-element energies to data
	std::vector< Triplet<double> > trips;
	append_energies(options,data,trips);
	int n_row_D = trips.back().row()+1;
	double dt2 = options->timestep_s * options->timestep_s;
	if (dt2 <= 0)
		dt2 = 1.0; // static solve

	// Weight matrix
	RowSparseMatrix<double> W2(n_row_D,n_row_D);
	VectorXi W_nnz = VectorXi::Ones(n_row_D);
	W2.reserve(W_nnz);
	int ne = data->indices.size();
	for (int i=0; i<ne; ++i)
	{
		const Vector2i &idx = data->indices[i];
		for (int j=0; j<idx[1]; ++j)
		{
			W2.coeffRef(idx[0]+j,idx[0]+j) = data->weights[i]*data->weights[i];
		}
	}

	// Weighted Laplacian
	data->D.resize(n_row_D,nx);
	data->D.setFromTriplets(trips.begin(), trips.end());
	data->Dt = data->D.transpose();
	data->DtW2 = dt2 * data->Dt * W2;
	data->A = data->DtW2 * data->D;
	for (int i=0; i<nx; ++i)
		data->A.coeffRef(i,i) += data->m[i];

	data->ldltA.compute(data->A);
	data->b.resize(nx,3);
	data->b.setZero();

	data->spring_k = options->mult_k*data->A.diagonal().maxCoeff();
	data->l = VectorXd::Zero(1);
	for (int i=0; i<3; ++i)
		data->K[i].resize(1,nx);

	// ADMM variables
	data->z.resize(n_row_D,3);
	data->z.setZero();
	data->u.resize(n_row_D,3);
	data->u.setZero();

} // end compute matrices

void Solver::append_energies(
	const Options *options,
	Data *data,
	std::vector<Triplet<double> > &D_triplets)
{
	int nt = data->tets.rows();
	if (nt==0)
		return;

	data->indices.reserve(nt);
	data->rest_volumes.reserve(nt);
	data->weights.reserve(nt);
	Lame lame;

	int energy_index = 0;
	for (int i=0; i<nt; ++i)
	{
		RowVector4i ele = data->tets.row(i);

		data->rest_volumes.emplace_back();
		data->weights.emplace_back();
		int energy_dim = EnergyTerm().init_tet(
			energy_index,
			lame,
			ele,
			&data->x,
			data->rest_volumes.back(),
			data->weights.back(),
			D_triplets );

		// Error in initialization
		if( energy_dim <= 0 ){
			data->rest_volumes.pop_back();
			data->weights.pop_back();
			continue;
		}

		data->indices.emplace_back(energy_index, energy_dim);
		energy_index += energy_dim;
	}
} // end append energies

} // namespace admmpd