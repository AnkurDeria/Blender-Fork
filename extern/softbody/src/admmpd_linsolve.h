// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_LINSOLVE_H_
#define ADMMPD_LINSOLVE_H_

#include "admmpd_types.h"

namespace admmpd {

// Preconditioned Conjugate Gradients
class ConjugateGradients {
public:
	void solve(
		const Options *options,
		SolverData *data);

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
		SolverData *data);

protected:
	void init_solve(
		const Options *options,
		SolverData *data);

	void compute_colors(
		const RowSparseMatrix<double> *A,
		int stride,
		std::vector<std::vector<int> > &colors);

};

} // namespace admmpd

#endif // ADMMPD_LINSOLVE_H_
