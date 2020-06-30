// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_LINSOLVE_H_
#define ADMMPD_LINSOLVE_H_

#include "admmpd_types.h"

namespace admmpd {

class GaussSeidel {
public:
	// Solves (A + KtK) x = (b + Ktl)
	// x and b passed as separate variables
	// for debugging/testing purposes.
	void solve(
		const Options *options,
		SolverData *data);

protected:
	// Allocates data, computes colors
	void init_solve(
		const Options *options,
		SolverData *data);

	// Computes colors of A + KtK
	void compute_colors(
		const RowSparseMatrix<double> *A,
		const RowSparseMatrix<double> *KtK, // if null, just A
		std::vector<std::vector<int> > &colors);

};

} // namespace admmpd

#endif // ADMMPD_LINSOLVE_H_
