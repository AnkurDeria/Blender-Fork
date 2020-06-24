// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_TETMESH_H_
#define ADMMPD_TETMESH_H_

#include "admmpd_types.h"

namespace admmpd {

class TetMesh {
public:
   // Given an embedding, compute masses
   // for the lattice vertices
   static void compute_masses(
        TetMeshData *mesh,
        const Eigen::MatrixXd *x,
        Eigen::VectorXd *masses,
        double density_kgm3 = 1100);

}; // class Lattice

} // namespace admmpd

#endif // ADMMPD_LATTICE_H_
