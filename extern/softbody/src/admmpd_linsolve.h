// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_LINSOLVE_H_
#define ADMMPD_LINSOLVE_H_

#include "admmpd_types.h"
#include "admmpd_collision.h"

namespace admmpd {

// Preconditioned Conjugate Gradients
class ConjugateGradients {
public:
	void solve(
		const Options *options,
		SolverData *data,
		Collision *collision);

protected:
	// Apply preconditioner
	void solve_Ax_b(
		SolverData *data,
		Eigen::VectorXd *x,
		Eigen::VectorXd *b);
};

// Multi-Colored Gauss-Seidel
class GaussSeidel {
public:
	void solve(
		const Options *options,
		SolverData *data,
		Collision *collision);

protected:
	void init_solve(
		const Options *options,
		SolverData *data,
		Collision *collision);

	void compute_colors(
		const std::vector<std::set<int> > &vertex_energies_graph,
		const std::vector<std::set<int> > &vertex_constraints_graph,
		std::vector<std::vector<int> > &colors);

	// For debugging:
	void verify_colors(SolverData *data);

};

} // namespace admmpd

#endif // ADMMPD_LINSOLVE_H_
