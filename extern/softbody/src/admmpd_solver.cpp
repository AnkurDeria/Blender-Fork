// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_solver.h"
#include "admmpd_energy.h"
#include "admmpd_collision.h"

#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <stdio.h>
#include <iostream>
#include <unordered_map>

#include "BLI_task.h" // threading
#include "BLI_assert.h"

namespace admmpd {
using namespace Eigen;
template <typename T> using RowSparseMatrix = SparseMatrix<T,RowMajor>;

typedef struct ThreadData {
	const Options *options;
	SolverData *data;
} ThreadData;

bool Solver::init(
    const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &T,
    const Options *options,
    SolverData *data)
{
	BLI_assert(data != NULL);
	BLI_assert(options != NULL);
	BLI_assert(V.rows() > 0);
	BLI_assert(V.cols() == 3);
	BLI_assert(T.rows() > 0);
	BLI_assert(T.cols() == 4);

	data->x = V;
	data->v.resize(V.rows(),3);
	data->v.setZero();
	data->tets = T;
	if (!compute_matrices(options,data))
		return false;

	return true;
} // end init

int Solver::solve(
	const Options *options,
	SolverData *data,
	Collision *collision)
{
	BLI_assert(data != NULL);
	BLI_assert(options != NULL);
	BLI_assert(data->x.cols() == 3);
	BLI_assert(data->x.rows() > 0);
	BLI_assert(data->A.nonZeros() > 0);
	BLI_assert(options->max_admm_iters > 0);

	// Init the solve which computes
	// quantaties like M_xbar and makes sure
	// the variables are sized correctly.
	init_solve(options,data);

	// Begin solver loop
	int iters = 0;
	for (; iters < options->max_admm_iters; ++iters)
	{
		// Update ADMM z/u
		solve_local_step(options,data);

		// Perform collision detection and linearization
		update_constraints(options,data,collision);

		// Solve Ax=b s.t. Kx=l
		data->b.noalias() = data->M_xbar + data->DtW2*(data->z-data->u);
		solve_conjugate_gradients(options,data);

	} // end solver iters

	// Update velocity (if not static solve)
	double dt = options->timestep_s;
	if (dt > 0.0)
		data->v.noalias() = (data->x-data->x_start)*(1.0/dt);

	return iters;
} // end solve

void Solver::init_solve(
	const Options *options,
	SolverData *data)
{
	BLI_assert(data != NULL);
	BLI_assert(options != NULL);
	int nx = data->x.rows();
	BLI_assert(nx > 0);

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

static void parallel_zu_update(
	void *__restrict userdata,
	const int i,
	const TaskParallelTLS *__restrict UNUSED(tls))
{
	ThreadData *td = (ThreadData*)userdata;
	Lame lame;
	lame.set_from_youngs_poisson(td->options->youngs,td->options->poisson);
	EnergyTerm().update(
		td->data->indices[i][0],
		lame,
		td->data->rest_volumes[i],
		td->data->weights[i],
		&td->data->x,
		&td->data->Dx,
		&td->data->z,
		&td->data->u );
} // end parallel zu update

void Solver::solve_local_step(
	const Options *options,
	SolverData *data)
{
	BLI_assert(data != NULL);
	BLI_assert(options != NULL);
	int ne = data->rest_volumes.size();
	BLI_assert(ne > 0);
  	ThreadData thread_data = {.options=options, .data = data};
	TaskParallelSettings settings;
	BLI_parallel_range_settings_defaults(&settings);
	BLI_task_parallel_range(0, ne, &thread_data, parallel_zu_update, &settings);
} // end local step

void Solver::update_constraints(
	const Options *options,
	SolverData *data,
	Collision *collision)
{
	BLI_assert(data != NULL);
	BLI_assert(options != NULL);

	if (collision==NULL)
		return;

	collision->detect(&data->x_start, &data->x);

	std::vector<double> l_coeffs;
	std::vector<Eigen::Triplet<double> > trips_x;
    std::vector<Eigen::Triplet<double> > trips_y;
    std::vector<Eigen::Triplet<double> > trips_z;

	// TODO collision detection
//	collision->jacobian(
//		&data->x,
//		&trips_x,
//		&trips_y,
//		&trips_z,
//		&l_coeffs);

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

typedef struct LinSolveThreadData {
	SolverData *data;
	MatrixXd *ls_x;
	MatrixXd *ls_b;
} LinSolveThreadData;

static void parallel_lin_solve(
	void *__restrict userdata,
	const int i,
	const TaskParallelTLS *__restrict UNUSED(tls))
{
	LinSolveThreadData *td = (LinSolveThreadData*)userdata;
	td->ls_x->col(i) = td->data->ldltA.solve(td->ls_b->col(i));
} // end parallel lin solve

void Solver::solve_conjugate_gradients(
	const Options *options,
	SolverData *data)
{
	BLI_assert(data != NULL);
	BLI_assert(options != NULL);
	int nx = data->x.rows();
	BLI_assert(nx > 0);
	BLI_assert(data->b.rows() == nx);
	BLI_assert(data->A.rows() == nx);
	BLI_assert(data->A.cols() == nx);
	BLI_assert(data->K[0].cols() == nx);
	BLI_assert(data->K[1].cols() == nx);
	BLI_assert(data->K[2].cols() == nx);
	BLI_assert(data->l.rows() > 0);
	BLI_assert(data->K[0].rows() == data->l.rows());
	BLI_assert(data->K[1].rows() == data->l.rows());
	BLI_assert(data->K[2].rows() == data->l.rows());

	// Solve Ax = b in parallel
	auto solve_Ax_b = [](
		SolverData *data_,
		MatrixXd *x_,
		MatrixXd *b_)
	{
		LinSolveThreadData thread_data = {.data=data_, .ls_x=x_, .ls_b=b_};
		TaskParallelSettings settings;
		BLI_parallel_range_settings_defaults(&settings);
		BLI_task_parallel_range(0, 3, &thread_data, parallel_lin_solve, &settings);
	};

	// If we don't have any constraints,
	// we don't need to perform CG
	if (std::max(std::max(
		data->K[0].nonZeros(),
		data->K[1].nonZeros()),
		data->K[2].nonZeros())==0)
	{
		solve_Ax_b(data,&data->x,&data->b);
		return;
	}

	// Inner product of matrices interpreted
	// if they were instead vectorized
	auto mat_inner = [](
		const MatrixXd &A,
		const MatrixXd &B)
	{
		double dot = 0.0;
		int nr = std::min(A.rows(), B.rows());
		for( int i=0; i<nr; ++i )
			for(int j=0; j<3; ++j)
				dot += A(i,j)*B(i,j);

		return dot;
	};

	// Update CGData
	admmpd::SolverData::CGData *cgdata = &data->cgdata;
	double eps = options->min_res;
	cgdata->b = data->b;
	if (cgdata->r.rows() != nx)
	{
		cgdata->r.resize(nx,3);
		cgdata->z.resize(nx,3);
		cgdata->p.resize(nx,3);
		cgdata->Ap.resize(nx,3);
	}

	for (int i=0; i<3; ++i)
	{
		RowSparseMatrix<double> Kt = data->K[i].transpose();
		cgdata->A[i] = data->A + data->spring_k*RowSparseMatrix<double>(Kt*data->K[i]);
		cgdata->b.col(i).noalias() += data->spring_k*Kt*data->l;
		cgdata->r.col(i).noalias() = cgdata->b.col(i) - cgdata->A[i]*data->x.col(i);
	}
	solve_Ax_b(data,&cgdata->z,&cgdata->r);
	cgdata->p = cgdata->z;

	for (int iter=0; iter<options->max_cg_iters; ++iter)
	{
		for( int i=0; i<3; ++i )
			cgdata->Ap.col(i).noalias() = cgdata->A[i]*cgdata->p.col(i);

		double p_dot_Ap = mat_inner(cgdata->p,cgdata->Ap);
		if( p_dot_Ap==0.0 )
			break;

		double zk_dot_rk = mat_inner(cgdata->z,cgdata->r);
		if( zk_dot_rk==0.0 )
			break;

		double alpha = zk_dot_rk / p_dot_Ap;
		data->x.noalias() += alpha * cgdata->p;
		cgdata->r.noalias() -= alpha * cgdata->Ap;
		if( cgdata->r.lpNorm<Infinity>() < eps )
			break;

		solve_Ax_b(data,&cgdata->z,&cgdata->r);
		double beta = mat_inner(cgdata->z,cgdata->r) / zk_dot_rk;
		cgdata->p = cgdata->z + beta*cgdata->p;
	}

} // end solve conjugate gradients

bool Solver::compute_matrices(
	const Options *options,
	SolverData *data)
{
	BLI_assert(data != NULL);
	BLI_assert(options != NULL);
	int nx = data->x.rows();
	BLI_assert(nx > 0);
	BLI_assert(data->x.cols() == 3);

	// Allocate per-vertex data
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
		compute_masses(options,data);

	// Add per-element energies to data
	std::vector<Triplet<double> > trips;
	append_energies(options,data,trips);
	if (trips.size()==0)
	{
		printf("**admmpd::Solver Error: No reduction coeffs\n");
		return false;
	}
	int n_row_D = trips.back().row()+1;
	double dt2 = options->timestep_s * options->timestep_s;
	if (options->timestep_s <= 0)
		dt2 = 1.0; // static solve, use dt=1 to not scale matrices

	// Diagonal weight matrix
	RowSparseMatrix<double> W2(n_row_D,n_row_D);
	VectorXi W_nnz = VectorXi::Ones(n_row_D);
	W2.reserve(W_nnz);
	int ne = data->indices.size();
	for (int i=0; i<ne; ++i)
	{
		const Vector2i &idx = data->indices[i];
		for (int j=0; j<idx[1]; ++j)
			W2.coeffRef(idx[0]+j,idx[0]+j) = data->weights[i]*data->weights[i];
	}

	// Mass weighted Laplacian
	data->D.resize(n_row_D,nx);
	data->D.setFromTriplets(trips.begin(), trips.end());
	data->DtW2 = dt2 * data->D.transpose() * W2;
	data->A = data->DtW2 * data->D;
	for (int i=0; i<nx; ++i)
		data->A.coeffRef(i,i) += data->m[i];
	data->ldltA.compute(data->A);
	data->b.resize(nx,3);
	data->b.setZero();

	// Constraint data
	data->spring_k = options->mult_k*data->A.diagonal().maxCoeff();
	data->l = VectorXd::Zero(1);
	for (int i=0; i<3; ++i)
		data->K[i].resize(1,nx);

	// ADMM dual/lagrange
	data->z.resize(n_row_D,3);
	data->z.setZero();
	data->u.resize(n_row_D,3);
	data->u.setZero();

	return true;

} // end compute matrices

void Solver::compute_masses(
	const Options *options,
	SolverData *data)
{
	BLI_assert(data != NULL);
	BLI_assert(options != NULL);

	// Source: https://github.com/mattoverby/mclscene/blob/master/include/MCL/TetMesh.hpp
	// Computes volume-weighted masses for each vertex
	// density_kgm3 is the unit-volume density
	data->m.resize(data->x.rows());
	data->m.setZero();
	int n_tets = data->tets.rows();
	for (int t=0; t<n_tets; ++t)
	{
		RowVector4i tet = data->tets.row(t);
		RowVector3d tet0 = data->x.row(tet[0]);
		Matrix3d edges;
		edges.col(0) = data->x.row(tet[1]) - tet0;
		edges.col(1) = data->x.row(tet[2]) - tet0;
		edges.col(2) = data->x.row(tet[3]) - tet0;
		double v = std::abs((edges).determinant()/6.f);
		double tet_mass = options->density_kgm3 * v;
		data->m[tet[0]] += tet_mass / 4.f;
		data->m[tet[1]] += tet_mass / 4.f;
		data->m[tet[2]] += tet_mass / 4.f;
		data->m[tet[3]] += tet_mass / 4.f;
	}
	// Verify masses
	int nx = data->m.rows();
	for (int i=0; i<nx; ++i)
	{
		if (data->m[i] <= 0.0)
		{
			printf("**Solver::compute_masses Error: unreferenced vertex\n");
			data->m[i]=1;
		}
	}
} // end compute masses

void Solver::append_energies(
	const Options *options,
	SolverData *data,
	std::vector<Triplet<double> > &D_triplets)
{
	BLI_assert(data != NULL);
	BLI_assert(options != NULL);
	int nt = data->tets.rows();
	BLI_assert(nt > 0);

	data->indices.reserve(nt);
	data->rest_volumes.reserve(nt);
	data->weights.reserve(nt);
	Lame lame;
	lame.set_from_youngs_poisson(options->youngs, options->poisson);

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
