


#ifndef _ADMMPD_H
#define _ADMMPD_H

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <vector>

namespace admmpd {

struct ADMMPD_Options {
    double timestep_s;
    int max_iters;
    Eigen::Vector3d grav;
    ADMMPD_Options() :
        timestep_s(1.0/100.0), // TODO: Figure out delta time from blender api!
        max_iters(20),
        grav(0,0,-9.8)
        {}
};

struct ADMMPD_Data {
    // Input:
    Eigen::MatrixXi tets; // elements lattice, t x 4
    Eigen::MatrixXd x; // vertices of lattice, n x 3
    // Set in compute_matrices: 
    Eigen::MatrixXd x_start; // x at beginning of timestep, n x 3
    Eigen::MatrixXd v; // velocity of lattice mesh, n x 3
    Eigen::VectorXd m; // masses of lattice verts, n x 1
    Eigen::MatrixXd z, z_prev; // ADMM z variable
    Eigen::MatrixXd u, u_prev; // ADMM u aug lag with W inv
    Eigen::MatrixXd M_xbar; // M*(x + dt v)
    Eigen::MatrixXd Dx; // D * x
    Eigen::MatrixXd b; // M xbar + DtW2(z-u)
    template <typename T> using RowSparseMatrix = Eigen::SparseMatrix<T,Eigen::RowMajor>;
    RowSparseMatrix<double> D; // reduction matrix
    RowSparseMatrix<double> Dt; // transpose reduction matrix
    RowSparseMatrix<double> DtW2; // D'W^2
    RowSparseMatrix<double> A; // M + D'W^2D
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ldltA;
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
        const ADMMPD_Options *options,
        ADMMPD_Data *data);

    // Solve a single time step.
    // Returns number of iterations.
    int solve(
        const ADMMPD_Options *options,
        ADMMPD_Data *data);

protected:

    void init_solve(
        const ADMMPD_Options *options,
        ADMMPD_Data *data);

    void compute_lattice(
        const ADMMPD_Options *options,
        ADMMPD_Data *data);

    void compute_matrices(
        const ADMMPD_Options *options,
        ADMMPD_Data *data);

	void append_energies(
		const ADMMPD_Options *options,
		ADMMPD_Data *data,
		std::vector<Eigen::Triplet<double> > &D_triplets);

}; // class ADMMPD_solver

} // namespace admmpd

#endif // _ADMMPD_H
