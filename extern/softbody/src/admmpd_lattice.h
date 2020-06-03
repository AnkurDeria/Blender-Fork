


#ifndef _ADMMPD_LATTICE_H
#define _ADMMPD_LATTICE_H

#include <Eigen/Dense>
#include <vector>

namespace admmpd {

class Lattice {
public:

    Eigen::MatrixXd vtx; // rest embedded vertices, p x 3
    Eigen::VectorXi vtx_to_tet; // what tet vtx is embedded in, p x 1
    Eigen::MatrixXd barys; // barycoords of the embedding, p x 4

    // Returns true on success
    bool generate(
        const Eigen::MatrixXd &V,
        Eigen::MatrixXd *x, // lattice vertices, n x 3
        Eigen::MatrixXi *tets); // lattice elements, m x 4

    // Given boxmin and boxmax, adds
    // 5 tets (and verts) to the vectors
    void create_packed_tets(
        const Eigen::Vector3d &min,
        const Eigen::Vector3d &max,
        std::vector<Eigen::Vector3d> &verts,
        std::vector<Eigen::Vector4i> &tets );

    // Returns true on success
    // Works on ptrs intead of local vars in case
    // I want to use it for something else.
    bool compute_vtx_tet_mapping(
        const Eigen::MatrixXd *vtx_, // embedded vertices, p x 3
        Eigen::VectorXi *vtx_to_tet_, // what tet vtx is embedded in, p x 1
        Eigen::MatrixXd *barys_, // barycoords of the embedding, p x 4
        const Eigen::MatrixXd *x_, // lattice vertices, n x 3
        const Eigen::MatrixXi *tets_); // lattice elements, m x 4

    // Returns the vtx mapped from x/v and tets
    Eigen::Vector3d get_mapped_vertex(
        int idx,
        const Eigen::MatrixXd *x_or_v,
        const Eigen::MatrixXi *tets );
//    void map_to_object(
//        const Eigen::MatrixXd *x,
//        const Eigen::MatrixXi *tets,
//        float (*vertexCos)[3]);

}; // class Lattice

} // namespace admmpd

#endif // _ADMMPD_LATTICE_H