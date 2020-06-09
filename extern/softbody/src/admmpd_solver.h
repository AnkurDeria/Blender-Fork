


#ifndef _ADMMPD_H
#define _ADMMPD_H

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <vector>

namespace admmpd {

struct Options {
    double timestep_s;
    int max_admm_iters;
    int max_cg_iters;
    double mult_k; // stiffness multiplier for constraints
    double min_res; // min residual for CG solver
    Eigen::Vector3d grav;
    Options() :
        timestep_s(1.0/100.0), // TODO: Figure out delta time from blender api
        max_admm_iters(20),
        max_cg_iters(10),
        mult_k(1.0),
        min_res(1e-4),
        grav(0,0,-9.8)
        {}
};

struct Data {
    // Input:
    Eigen::MatrixXi tets; // elements t x 4
    Eigen::MatrixXd x; // vertices, n x 3
    Eigen::MatrixXd v; // velocity, n x 3 TODO: from cache
    // Set in compute_matrices: 
    Eigen::MatrixXd x_start; // x at beginning of timestep, n x 3
    Eigen::VectorXd m; // masses, n x 1 TODO: from BodyPoint
    Eigen::MatrixXd z; // ADMM z variable
    Eigen::MatrixXd u; // ADMM u aug lag with W inv
    Eigen::MatrixXd M_xbar; // M*(x + dt v)
    Eigen::MatrixXd Dx; // D * x
    Eigen::MatrixXd b; // M xbar + DtW2(z-u)
    template <typename T> using RowSparseMatrix = Eigen::SparseMatrix<T,Eigen::RowMajor>;
    RowSparseMatrix<double> D; // reduction matrix
    RowSparseMatrix<double> Dt; // transpose reduction matrix
    RowSparseMatrix<double> DtW2; // D'W^2
    RowSparseMatrix<double> A; // M + D'W^2D
    RowSparseMatrix<double> K[3]; // constraint Jacobian
    Eigen::VectorXd l; // constraint rhs (Kx=l)
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ldltA;
    double spring_k;
    // Set in append_energies:
	std::vector<Eigen::Vector2i> indices; // per-energy index into D (row, num rows)
	std::vector<double> rest_volumes; // per-energy rest volume
	std::vector<double> weights; // per-energy weights
};

class Solver {
public:

    // Initialies solver data. If a per-vertex
    // variable is resized it is initialized to zero.
    // Returns true on success
    bool init(
        const Eigen::MatrixXd &V, // vertices
        const Eigen::MatrixXi &T, // tets
        const Options *options,
        Data *data);

    // Solve a single time step.
    // Returns number of iterations.
    int solve(
        const Options *options,
        Data *data);

protected:

    void update_constraints(
        const Options *options,
        Data *data);

    void init_solve(
        const Options *options,
        Data *data);

    // Global step with CG:
    // 1/2||Ax-b||^2 + k/2||Kx-l||^2
	void solve_conjugate_gradients(
        const Options *options,
        Data *data);

    void compute_lattice(
        const Options *options,
        Data *data);

    void compute_matrices(
        const Options *options,
        Data *data);

	void append_energies(
		const Options *options,
		Data *data,
		std::vector<Eigen::Triplet<double> > &D_triplets);

}; // class ADMMPD_solver

} // namespace admmpd

#endif // _ADMMPD_H
