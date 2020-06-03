


#include "admmpd_solver.h"
#include "admmpd_lattice.h"
#include "admmpd_energy.h"

#include <Eigen/Geometry>
#include <Eigen/Sparse>

#include <stdio.h>

namespace admmpd {
using namespace Eigen;
template <typename T> using RowSparseMatrix = SparseMatrix<T,RowMajor>;

bool Solver::init(
    const Eigen::MatrixXd &V,
	const Eigen::MatrixXi &T,
    const ADMMPD_Options *options,
    ADMMPD_Data *data)
{
	if (!data || !options)
		throw std::runtime_error("init: data/options null");

	data->x = V;
	data->tets = T;
	compute_matrices(options,data);
	return true;
} // end init

void Solver::init_solve(
	const ADMMPD_Options *options,
	ADMMPD_Data *data)
{
	double dt = std::max(0.0, options->timestep_s);
	int nx = data->x.rows();
	if (data->M_xbar.rows() != nx)
		data->M_xbar.resize(nx,3);

	// velocity and position
	data->x_start = data->x;
	for( int i=0; i<nx; ++i )
	{
		data->v.row(i) += options->grav;
		data->M_xbar.row(i) =
			data->m[i] * data->x.row(i) +
			dt*data->m[i]*data->v.row(i);
	}

	// ADMM variables
	data->Dx.noalias() = data->D * data->x;
	data->z = data->Dx;
	data->z_prev = data->z;
	data->u.setZero();
	data->u_prev.setZero();

} // end init solve

int Solver::solve(
	const ADMMPD_Options *options,
	ADMMPD_Data *data)
{
	// Init the solve which computes
	// quantaties like M_xbar and makes sure
	// the variables are sized correctly.
	init_solve(options,data);

	int ne = data->rest_volumes.size();
	Lame lame;

	// Begin solver loop
	int iters = 0;
	int max_iters = options->max_iters < 0 ? 10000 : options->max_iters;
	for (; iters < max_iters; ++iters)
	{
		// Local step
		data->Dx.noalias() = data->D * data->x;
		data->z_prev.noalias() = data->z;
		data->u_prev.noalias() = data->u;
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
		data->b.noalias() = data->M_xbar + data->DtW2*(data->z-data->u);
		data->x.noalias() = data->ldltA.solve(data->b);

	} // end solver iters

	double dt = std::max(0.0, options->timestep_s);
	if (dt > 0.0)
		data->v.noalias() = (data->x-data->x_start)*(1.0/dt);

	return iters;
} // end solve

void Solver::compute_matrices(
	const ADMMPD_Options *options,
	ADMMPD_Data *data)
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
	}

	// Add per-element energies to data
	std::vector< Triplet<double> > trips;
	append_energies(options,data,trips);
	int nw = trips.back().row()+1;

	// Global matrix
	data->D.resize(nw,nx);
	data->D.setFromTriplets(trips.begin(), trips.end());
	data->Dt = data->D.transpose();
	VectorXd w_diag = Map<VectorXd>(data->weights.data(), data->weights.size());
	data->DtW2 = data->Dt * w_diag.asDiagonal() * w_diag.asDiagonal();
	data->A = data->DtW2 * data->D;
	data->A.diagonal() += data->m;
	data->ldltA.compute(data->A);
	data->b.resize(nx,3);
	data->b.setZero();

	data->z.resize(nw,3);
	data->z.setZero();
	data->z_prev.resize(nw,3);
	data->z_prev.setZero();
	data->u.resize(nw,3);
	data->u.setZero();
	data->u_prev.resize(nw,3);
	data->u_prev.setZero();

} // end compute matrices

void Solver::append_energies(
	const ADMMPD_Options *options,
	ADMMPD_Data *data,
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
		data->rest_volumes.emplace_back();
		data->weights.emplace_back();
		int energy_dim = 0;

		RowVector4i ele = data->tets.row(i);
		energy_dim = EnergyTerm().init_tet(
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
