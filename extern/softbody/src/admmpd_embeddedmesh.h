// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_EMBEDDEDMESH_H_
#define ADMMPD_EMBEDDEDMESH_H_

#include "admmpd_types.h"

namespace admmpd {

class EmbeddedMesh {
public:
    // Returns true on success
    bool generate(
        const Eigen::MatrixXd &V, // embedded verts
        const Eigen::MatrixXi &F, // embedded faces
        EmbeddedMeshData *emb_mesh, // where embedding is stored
        Eigen::MatrixXd *x_tets); // lattice vertices, n x 3

    // Returns the vtx mapped from x/v and tets
    Eigen::Vector3d get_mapped_vertex(
        const EmbeddedMeshData *emb_mesh,
        const Eigen::MatrixXd *x_data,
        int idx);

    // Given an embedding, compute masses
    // for the lattice vertices
    void compute_masses(
        EmbeddedMeshData *emb_mesh, // where embedding is stored
        const Eigen::MatrixXd *x_embed, // embedded vertices, p x 3
        const Eigen::MatrixXd *x_tets, // lattice vertices, n x 3
        Eigen::VectorXd *masses_tets, // masses of the lattice verts
        double density_kgm3 = 1100);

protected:

    // Returns true on success
    // Computes the embedding data, like barycoords
    bool compute_embedding(
        EmbeddedMeshData *emb_mesh, // where embedding is stored
        const Eigen::MatrixXd *x_embed, // embedded vertices, p x 3
        const Eigen::MatrixXd *x_tets); // lattice vertices, n x 3

}; // class Lattice

} // namespace admmpd

#endif // ADMMPD_EMBEDDEDMESH_H_
