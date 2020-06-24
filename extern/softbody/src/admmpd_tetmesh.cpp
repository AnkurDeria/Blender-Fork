// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_tetmesh.h"
#include "admmpd_math.h"
#include <unordered_map>
#include <set>
#include "BLI_assert.h"

namespace admmpd {
using namespace Eigen;

void TetMesh::compute_masses(
    TetMeshData *mesh,
    const Eigen::MatrixXd *x,
    Eigen::VectorXd *masses,
	double density_kgm3)
{
	BLI_assert(mesh != NULL);
	BLI_assert(x != NULL);
	BLI_assert(masses != NULL);
	BLI_assert(density_kgm3 > 0);

	// Source: https://github.com/mattoverby/mclscene/blob/master/include/MCL/TetMesh.hpp
	// Computes volume-weighted masses for each vertex
	// density_kgm3 is the unit-volume density
	int nx = x->rows();
	masses->resize(nx);
	masses->setZero();
	int n_tets = mesh->tets.rows();
	for (int t=0; t<n_tets; ++t)
	{
		RowVector4i tet = mesh->tets.row(t);
		RowVector3d tet_v0 = x->row(tet[0]);
		Matrix3d edges;
		edges.col(0) = x->row(tet[1]) - tet_v0;
		edges.col(1) = x->row(tet[2]) - tet_v0;
		edges.col(2) = x->row(tet[3]) - tet_v0;
		double vol = std::abs((edges).determinant()/6.f);
		double tet_mass = density_kgm3 * vol;
		masses->operator[](tet[0]) += tet_mass / 4.f;
		masses->operator[](tet[1]) += tet_mass / 4.f;
		masses->operator[](tet[2]) += tet_mass / 4.f;
		masses->operator[](tet[3]) += tet_mass / 4.f;
	}

	// Verify masses
	for (int i=0; i<nx; ++i)
	{
		if (masses->operator[](i) <= 0.0)
		{
			printf("**TetMesh::compute_masses Error: unreferenced vertex\n");
			masses->operator[](i)=1;
		}
	}
} // end compute masses

} // namespace admmpd
