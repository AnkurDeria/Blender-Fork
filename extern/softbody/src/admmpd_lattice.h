// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_LATTICE_H_
#define ADMMPD_LATTICE_H_

#include <Eigen/Dense>
#include <vector>

namespace admmpd {

class Lattice {
public:
    Eigen::VectorXi vtx_to_tet; // what tet vtx is embedded in, p x 1
    Eigen::MatrixXd barys; // barycoords of the embedding, p x 4

    // Returns true on success
    bool generate(
        const Eigen::MatrixXd &V,
        Eigen::MatrixXd *x, // lattice vertices, n x 3
        Eigen::MatrixXi *tets); // lattice elements, m x 4

    // Returns the vtx mapped from x/v and tets
    Eigen::Vector3d get_mapped_vertex(
        int idx,
        const Eigen::MatrixXd *x_or_v,
        const Eigen::MatrixXi *tets );

protected:

    // Returns true on success
    bool compute_vtx_tet_mapping(
        const Eigen::MatrixXd *vtx_, // embedded vertices, p x 3
        Eigen::VectorXi *vtx_to_tet_, // what tet vtx is embedded in, p x 1
        Eigen::MatrixXd *barys_, // barycoords of the embedding, p x 4
        const Eigen::MatrixXd *x_, // lattice vertices, n x 3
        const Eigen::MatrixXi *tets_); // lattice elements, m x 4

}; // class Lattice

} // namespace admmpd

#endif // ADMMPD_LATTICE_H_
