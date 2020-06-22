// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_collision.h"
#include "BLI_assert.h"

namespace admmpd {
using namespace Eigen;

void EmbeddedMeshCollision::set_obstacles(
	const float *v0,
	const float *v1,
	int nv,
	const int *faces,
	int nf)
{
	if (obs_V0.rows() != nv)
		obs_V0.resize(nv,3);

	if (obs_V1.rows() != nv)
		obs_V1.resize(nv,3);

	for (int i=0; i<nv; ++i)
	{
		for (int j=0; j<3; ++j)
		{
			obs_V0(i,j) = v0[i*3+j];
			obs_V1(i,j) = v1[i*3+j];
		}
	}

	if (obs_F.rows() != nf)
	{
		obs_F.resize(nf,3);
		obs_aabbs.resize(nf);
	}

	for (int i=0; i<nf; ++i)
	{
		obs_aabbs[i].setEmpty();
		for (int j=0; j<3; ++j)
		{
			int fj = faces[i*3+j];
			obs_F(i,j) = fj;
			obs_aabbs[i].extend(obs_V0.row(fj).transpose());
			obs_aabbs[i].extend(obs_V1.row(fj).transpose());
		}
	}

	obs_tree.init(obs_aabbs);

} // end add obstacle

void EmbeddedMeshCollision::detect(const Eigen::MatrixXd *x0, const Eigen::MatrixXd *x1)
{
	/*
	// First, update the positions of the embedded vertex
	// and perform collision detection against obstacles
	int n_ev = emb_V0.rows();
	std::vector<int> colliding(n_ev,0);
	for (int i=0; i<n_ev; ++i)
	{
		int t = vtx_to_tet->operator[](i);
		BLI_assert(t >= 0);
		BLI_assert(t < tets->rows());
		RowVector4i tet = tets->row(t);
		Vector4d bary = emb_barys->row(i);
//		emb_V0.row(i) =
//			bary[0] * x0->row(tet[0]) +
//			bary[1] * x0->row(tet[1]) +
//			bary[2] * x0->row(tet[2]) +
//			bary[3] * x0->row(tet[3]);
		Vector3d pt = 
			bary[0] * x1->row(tet[0]) +
			bary[1] * x1->row(tet[1]) +
			bary[2] * x1->row(tet[2]) +
			bary[3] * x1->row(tet[3]);
//		emb_V1.row(i) =

		// Check if we are inside the mesh.
		// If so, find the nearest face in the rest pose.
		PointInTriangleMeshTraverse<double> pt_in_mesh(pt, &obs_V1, &obs_F);
		obs_tree.traverse(pt_in_mesh);
		if (pt_in_mesh.output.num_hits() % 2 == 1)
		{
			// Find 
//			hit = true;
		}

//		colliding[i] = hit;
	}

	// Only bother with self collision if it
	// is not colliding with an obstacle.
	// This is only useful for discrete tests.
*/
} // end emb collision detect

/*
void FloorCollider::detect(const Eigen::MatrixXd *x)
{
	(void)(x);
	// Can just do detection in jacobian I guess.
}

void FloorCollider::jacobian(
	const Eigen::MatrixXd *x,
	std::vector<Eigen::Triplet<double> > *trips_x,
	std::vector<Eigen::Triplet<double> > *trips_y,
	std::vector<Eigen::Triplet<double> > *trips_z,
	std::vector<double> *l)
{
	const double floor_z = -2.0;

	int nx = x->rows();
	for (int i=0; i<nx; ++i)
	{
		Eigen::Vector3d xi = x->row(i);
		if (xi[2]>floor_z)
			continue;

		int c_idx = l->size();
		trips_z->emplace_back(c_idx,i,1.0);
		l->emplace_back(floor_z);
	}
} // end floor collider Jacobian
*/
} // namespace admmpd
