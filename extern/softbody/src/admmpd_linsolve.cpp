// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_linsolve.h"
#include <numeric>
#include <iostream>
#include <unordered_set>
#include "BLI_assert.h"

namespace admmpd {
using namespace Eigen;

void GaussSeidel::solve(
		const Options *options,
		SolverData *data)
{
	init_solve(options,data);
	std::vector<std::vector<int> > *colors;
	//if (data->gsdata.KtK.nonZeros()==0)
		colors = &data->gsdata.A_colors;
	//else...

	double omega = 1.0; // over relaxation
	int n_colors = colors->size();

	// Outer iteration loop
	int iter = 0;
	for (; iter < options->max_gs_iters; ++iter)
	{
		for (int color=0; color<n_colors; ++color)
		{
			const std::vector<int> &inds = colors->at(color);
			int n_inds = inds.size();
			for (int i=0; i<n_inds; ++i)
			{
				int idx = inds[i];

				// Special case pins TODO
				// We can skip the usual Gauss-Seidel update
	//			if (is_pinned[idx]) ...

				RowSparseMatrix<double>::InnerIterator rit(data->A,idx);
				Vector3d LUx(0,0,0);
				Vector3d inv_aii(0,0,0);
				for (; rit; ++rit)
				{
					int r = rit.row();
					int c = rit.col();
					double v = rit.value();
					if (v==0.0)
						continue;

					if (r==c) // Diagonal
					{
						inv_aii.array() = 1.0/v;
						continue;
					}
					Vector3d xj = data->x.row(c);
					LUx += v*xj;
				}

				// Update x
				Vector3d bi = data->b.row(idx);
				Vector3d xi = data->x.row(idx);
				Vector3d xi_new = (bi-LUx);

				for (int j=0; j<3; ++j)
					xi_new[j] *= inv_aii[j];
				data->x.row(idx) = xi*(1.0-omega) + xi_new*omega;
				data->gsdata.last_dx.row(idx) = data->x.row(idx)-xi.transpose();

				// TODO
				// We can also apply constraints here, like
				// checking against Collision::floor_z
				if (data->x(idx,2)<0)
					data->x(idx,2)=0;

			} // end loop inds
		} // end loop colors

		// TODO check exit condition

	} // end loop GS iters

} // end solve with constraints

void GaussSeidel::init_solve(
		const Options *options,
		SolverData *data)
{
	BLI_assert(options != nullptr);
	BLI_assert(data != nullptr);
	int nx = data->x.rows();
	BLI_assert(nx>0);
	BLI_assert(data->x.cols()==3);
	data->b.noalias() = data->M_xbar + data->DtW2*(data->z-data->u);
	BLI_assert(data->b.rows()==nx);
	BLI_assert(data->b.cols()==data->x.cols());

	if (data->gsdata.last_dx.rows() != nx)
	{
		data->gsdata.last_dx.resize(nx,3);
		data->gsdata.last_dx.setZero();
	}

	// Do we need to color the default colorings?
	if( data->gsdata.A_colors.size() == 0 )
		compute_colors(&data->A,nullptr,data->gsdata.A_colors);

	// TODO: Eventually we'll replace KtK with the full-dof matrix.
	// For now use z and test collisions against ground plane.
	bool has_constraints = data->K[2].nonZeros()>0;
	data->gsdata.KtK = data->K[2].transpose()*data->K[2];

	if (false)//(has_constraints)
	{
		(void)(has_constraints);
		// TODO color A + KtK
	}

} // end init solve

void GaussSeidel::compute_colors(
		const RowSparseMatrix<double> *A,
		const RowSparseMatrix<double> *KtK, // if null, just A
		std::vector<std::vector<int> > &colors)
{
	BLI_assert(A != nullptr);
	if (KtK != nullptr)
	{
		throw std::runtime_error("**GaussSeidel::compute_colors TODO: KtK");
	}

	colors.clear();
	int nx = A->rows();

	{ // DEBUGGING
		colors.resize(nx,std::vector<int>());
		for (int i=0; i<nx; ++i)
		{
			colors[i].emplace_back(i);
		}
	}

} // end compute colors

} // namespace admmpd
