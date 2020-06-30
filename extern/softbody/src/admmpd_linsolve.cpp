// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_linsolve.h"
#include <numeric>
#include <iostream>
#include <set>
#include <unordered_set>
#include <random>
#include "BLI_assert.h"
#include "BLI_task.h"

namespace admmpd {
using namespace Eigen;

// Makes full n3 x n3 matrix
static inline void make_n3(
		const RowSparseMatrix<double> &A,
		RowSparseMatrix<double> &A3)
{
	int na = A.rows();
	SparseMatrix<double> A_cm = A; // A to CM
	SparseMatrix<double> A3_cm(na*3, na*3);
	A3_cm.reserve(A_cm.nonZeros()*3);
	int cols = A_cm.rows();
	for(int c=0; c<cols; ++c)
	{
		for(int dim=0; dim<3; ++dim)
		{
			int col = c*3+dim;
			A3_cm.startVec(col);
			for(SparseMatrix<double>::InnerIterator itA(A_cm,c); itA; ++itA)
				A3_cm.insertBack((itA.row()*3+dim), col) = itA.value();
		}
	}
	A3_cm.finalize();
	A3_cm.makeCompressed();
	A3 = A3_cm; // to RowMajor
} // end make_n3

void ConjugateGradients::solve(
		const Options *options,
		SolverData *data)
{
	BLI_assert(data != NULL);
	BLI_assert(options != NULL);
	int nx = data->x.rows();
	BLI_assert(nx > 0);
	BLI_assert(data->A.rows() == nx);
	BLI_assert(data->A.cols() == nx);
	BLI_assert(data->b.rows() == nx);
	BLI_assert(data->C.cols() == nx*3);
	BLI_assert(data->d.rows() > 0);
	BLI_assert(data->C.rows() == data->d.rows());

	// If we don't have any constraints, we don't need to perform CG
	data->b.noalias() = data->M_xbar + data->DtW2*(data->z-data->u);
	if (data->C.nonZeros()==0)
	{
		data->x = data->ldltA.solve(data->b);
		return;
	}

	// Resize data associated with Conjugate-Gradient
	admmpd::SolverData::GlobalStepData *gsdata = &data->gsdata;
	if (gsdata->r.rows() != nx*3)
	{
		gsdata->r.resize(nx*3);
		gsdata->z.resize(nx*3);
		gsdata->p.resize(nx*3);
		gsdata->A3p.resize(nx*3);
		gsdata->b3_plus_Ctd.resize(nx*3);
		make_n3(data->A, gsdata->A3);
	}

	gsdata->Ctd = data->spring_k * data->C.transpose()*data->d;
	gsdata->CtC = data->spring_k * data->C.transpose()*data->C;
	gsdata->A3_plus_CtC = gsdata->A3 + gsdata->CtC;
	VectorXd x3(nx*3);
	for (int i=0; i<nx; ++i)
	{
		Vector3d bi = data->b.row(i);
		gsdata->b3_plus_Ctd.segment<3>(i*3) = bi+gsdata->Ctd.segment<3>(i*3);
		x3.segment<3>(i*3) = data->x.row(i);
	}

	gsdata->r.noalias() = gsdata->b3_plus_Ctd - gsdata->A3_plus_CtC*x3;
	solve_Ax_b(data,&gsdata->z,&gsdata->r);
	gsdata->p = gsdata->z;

	for (int iter=0; iter<options->max_cg_iters; ++iter)
	{
		gsdata->A3p.noalias() = gsdata->A3_plus_CtC*gsdata->p;
		double p_dot_Ap = gsdata->p.dot(gsdata->A3p);
		if (p_dot_Ap==0.0)
			break;

		double zk_dot_rk = gsdata->z.dot(gsdata->r);
		if (zk_dot_rk==0.0)
			break;

		double alpha = zk_dot_rk / p_dot_Ap;
		x3.noalias() += alpha * gsdata->p;
		gsdata->r.noalias() -= alpha * gsdata->A3p;
		double r_norm = gsdata->r.lpNorm<Infinity>();
		if (r_norm < options->min_res)
			break;

		solve_Ax_b(data,&gsdata->z,&gsdata->r);
		double zk1_dot_rk1 = gsdata->z.dot(gsdata->r);
		double beta = zk1_dot_rk1 / zk_dot_rk;
		gsdata->p = gsdata->z + beta*gsdata->p;
	}

	for (int i=0; i<nx; ++i)
		data->x.row(i) = x3.segment<3>(i*3);

} // end ConjugateGradients solve

// Solve Ax = b in parallel and apply
// dimensionality mapping with temporaries
void ConjugateGradients::solve_Ax_b(
	SolverData *data_,
	VectorXd *x_,
	VectorXd *b_)
{
	typedef struct LinSolveThreadData {
		SolverData *data;
		VectorXd *ls_x;
		VectorXd *ls_b;
	} LinSolveThreadData;

	auto parallel_lin_solve = [](
		void *__restrict userdata,
		const int i,
		const TaskParallelTLS *__restrict UNUSED(tls))->void
	{
		LinSolveThreadData *td = (LinSolveThreadData*)userdata;
		int nx = td->ls_x->rows()/3;
		VectorXd b(nx);
		for (int j=0; j<nx; ++j)
			b[j] = td->ls_b->operator[](j*3+i);
		VectorXd x = td->data->ldltA.solve(b);
		for (int j=0; j<nx; ++j)
			td->ls_x->operator[](j*3+i) = x[j];
	};

	LinSolveThreadData thread_data = {.data=data_, .ls_x=x_, .ls_b=b_};
	TaskParallelSettings settings;
	BLI_parallel_range_settings_defaults(&settings);
	BLI_task_parallel_range(0, 3, &thread_data, parallel_lin_solve, &settings);

} // end solve Ax=b

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

	// Do we need to color the default colorings?
	if (data->gsdata.A_colors.size() == 0)
		compute_colors(&data->A, 1, data->gsdata.A_colors);

	// Verify color groups are correct
	{
		std::cout << "TESTING " << data->tets.rows() << " tets" << std::endl;
		std::cout << "num colors: " << data->gsdata.A_colors.size() << std::endl;
		int nt = data->tets.rows();
		int nc = data->gsdata.A_colors.size();
		for (int i=0; i<nc; ++i)
		{
			// Each vertex in the color should not
			// be a part of a tet with a vertex in the same color
			const std::vector<int> &grp = data->gsdata.A_colors[i];
			int n_grp = grp.size();
			for (int j=0; j<n_grp; ++j)
			{
				int p_idx = grp[j];
				auto in_tet = [](int idx, const RowVector4i &t)
				{
					return (t[0]==idx||t[1]==idx||t[2]==idx||t[3]==idx);
				};

				for (int k=0; k<nt; ++k)
				{
					RowVector4i t = data->tets.row(k);
					if (!in_tet(p_idx,t))
						continue;

					for (int l=0; l<n_grp; ++l)
					{
						int q_idx = grp[l];
						if (p_idx==q_idx)
							continue;

						if (in_tet(q_idx,t))
						{
std::cerr << "p: " << p_idx << ", q: " << q_idx << ", tet (" << k << "): " << t << std::endl;
//							throw std::runtime_error("GaussSeidel Error: Color contains two verts form same tet!");
						}
					}
				}
			}
		}
	}

	// Create large A if we haven't already.
	if (data->gsdata.A3.rows() != nx*3)
		make_n3(data->A, data->gsdata.A3);

	// TODO
	// Eventually we'll replace KtK with the full-dof matrix.
	// For now use z and test collisions against ground plane.
	bool has_constraints = data->C.nonZeros()>0;

	// Finally, the new global matrix and rhs
	if (has_constraints)
	{
		data->gsdata.CtC = data->spring_k * data->C.transpose()*data->C;
		data->gsdata.Ctd.noalias() = data->spring_k * data->C.transpose()*data->d;
		data->gsdata.A3_plus_CtC = data->gsdata.A3 + data->gsdata.CtC;
		data->gsdata.b3_plus_Ctd.resize(nx*3);
		for (int i=0; i<nx; ++i)
		{
			data->gsdata.b3_plus_Ctd[i*3+0] = data->b(i,0)+data->gsdata.Ctd[i*3+0];
			data->gsdata.b3_plus_Ctd[i*3+1] = data->b(i,1)+data->gsdata.Ctd[i*3+1];
			data->gsdata.b3_plus_Ctd[i*3+2] = data->b(i,2)+data->gsdata.Ctd[i*3+2];
		}
		compute_colors(&data->gsdata.A3_plus_CtC, 3, data->gsdata.A3_plus_CtC_colors);
	}
	else
	{
		data->gsdata.CtC.setZero();
		data->gsdata.Ctd.setZero();
		data->gsdata.A3_plus_CtC = data->gsdata.A3;
		data->gsdata.b3_plus_Ctd.resize(nx*3);
		for (int i=0; i<nx; ++i)
		{
			data->gsdata.b3_plus_Ctd[i*3+0] = data->b(i,0);
			data->gsdata.b3_plus_Ctd[i*3+1] = data->b(i,1);
			data->gsdata.b3_plus_Ctd[i*3+2] = data->b(i,2);
		}
		data->gsdata.A3_plus_CtC_colors = data->gsdata.A_colors;
	}

} // end init solve

typedef struct GraphColorThreadData {
	const RowSparseMatrix<double> *A;
	int stride;
	std::vector<std::vector<int> > *adjacency;
	std::vector<std::set<int> > *palette;
	int init_palette_size;
	std::vector<int> *conflict;
	std::vector<int> *node_colors;
} GraphColorThreadData;

// Rehash of graph coloring from
// https://github.com/mattoverby/mclscene/blob/master/include/MCL/GraphColor.hpp
void GaussSeidel::compute_colors(
		const RowSparseMatrix<double> *A,
		int stride,
		std::vector<std::vector<int> > &colors)
{
	BLI_assert(A != nullptr);
	BLI_assert(stride>0);

	// Graph color settings
	int init_palette_size = 6;

	// Graph coloring tmp data
	int n_nodes = A->rows()/stride;
	std::vector<std::vector<int> > adjacency(n_nodes, std::vector<int>());
	std::vector<std::set<int> > palette(n_nodes, std::set<int>());
	std::vector<int> conflict(n_nodes,1);
	std::vector<int> node_colors(n_nodes,-1);

	//
	// Step 1)
	// Graph color initialization
	//
	auto parallel_init_graph = [](
		void *__restrict userdata,
		const int i,
		const TaskParallelTLS *__restrict UNUSED(tls)) -> void
	{
		GraphColorThreadData *td = (GraphColorThreadData*)userdata;
		typedef RowSparseMatrix<double>::InnerIterator InnerIter;
		// Use lower trianglular portion of the matrix.
		// That is, store adjacency inds that are lower than
		// the current index.
		for (InnerIter it(*(td->A),i*td->stride); it; ++it)
		{
			if (std::abs(it.value())==0.0)
				continue;
			if (it.col()/td->stride >= it.row()/td->stride)
				break;
			int idx = it.col()/td->stride;
			td->adjacency->at(i).emplace_back(idx);
		}
		for( int j=0; j<td->init_palette_size; ++j ) // init colors
			td->palette->at(i).insert(j); 
	};
	GraphColorThreadData thread_data = {
		.A = A, .stride = stride,
		.adjacency = &adjacency,
		.palette = &palette,
		.init_palette_size = init_palette_size,
		.conflict = &conflict,
		.node_colors = &node_colors };
	TaskParallelSettings settings;
	BLI_parallel_range_settings_defaults(&settings);
	BLI_task_parallel_range(0, n_nodes, &thread_data, parallel_init_graph, &settings);

    std::random_device rd;
    std::mt19937 mt(rd());
	std::uniform_int_distribution<int> dist(0, n_nodes);

	//
	// Step 2)
	// Stochastic Graph coloring
	//
	std::vector<int> node_queue(n_nodes);
	std::iota(node_queue.begin(), node_queue.end(), 0);
	int max_iters = n_nodes;
	for (int rand_iter=0; n_nodes>0 && rand_iter < max_iters; ++rand_iter)
	{
		// Generate a random color and find conflicts
		for (int i=0; i<n_nodes; ++i)
		{
			int idx = node_queue[i];
			if( palette[idx].size() < 2 ){ // Feed if hungry
				palette[idx].insert(init_palette_size+rand_iter);
			}

			int c_idx = dist(mt) % (int)palette[idx].size();
			int curr_color = *std::next(palette[idx].begin(), c_idx);
			node_colors[idx] = curr_color;
			conflict[idx] = 0;

			// Hungarian heuristic: node with largest index keeps color
			// if both have the same color.
			// Check lower-indexed adjacent nodes, and mark them
			// in conflict if they generated the same color.
			// We can do this if we process colors in an increasing order.
			int n_adj = adjacency[idx].size();
			for (int j=0; j<n_adj; ++j)
			{
				// Adjacency index buffer only stores lower-indexed nodes.
				int adj_idx = adjacency[idx][j];
				int adj_color = node_colors[adj_idx];
				if (adj_idx >= idx)
					throw std::runtime_error("GaussSeidel Error: Oops, not lower diag");
				if (adj_color==curr_color)
					conflict[adj_idx] = 1;					
			}
		}

		// Resolve conflicts and trim queue
		std::set<int> new_queue;
		for (int i=0; i<n_nodes; ++i)
		{
			int idx = node_queue[i];
			int curr_color = node_colors[idx];
			int n_adj = adjacency[idx].size();

			// If this node not in conflict, remove its
			// color from adjacent palettes. If adjacent
			// nodes not in conflict, remove their color
			// from this palette.
			for (int j=0; j<n_adj; ++j)
			{
				int adj_idx = adjacency[idx][j];
				int adj_color = node_colors[adj_idx];
				if (!conflict[idx])
					palette[adj_idx].erase(curr_color);
				if (!conflict[adj_idx])
					palette[idx].erase(adj_color);
			}

			if (conflict[idx])
				new_queue.emplace(idx);
		}

		n_nodes = new_queue.size();
		node_queue.assign(new_queue.begin(), new_queue.end());

	} // end color loop

	//
	// Step 3)
	// Map per-vertex colors
	//
	colors.clear();
	n_nodes = node_colors.size();
	for( int i=0; i<n_nodes; ++i ){
		int color = node_colors[i];
		if (color<0)
			throw std::runtime_error("GaussSeidel: Error with coloring");
		while( color >= (int)colors.size() )
			colors.emplace_back(std::vector<int>());
		colors[color].emplace_back(i);
	}

	// Remove empty color groups
	for (std::vector<std::vector<int> >::iterator it = colors.begin(); it != colors.end();)
		it->size() == 0 ? it = colors.erase(it) : it++;

} // end compute colors

} // namespace admmpd