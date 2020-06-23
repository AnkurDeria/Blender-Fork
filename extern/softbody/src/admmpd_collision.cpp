// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_collision.h"
#include "admmpd_bvh_traverse.h"
#include "BLI_assert.h"
#include "BLI_task.h"
#include "BLI_threads.h"

namespace admmpd {
using namespace Eigen;

VFCollisionPair::VFCollisionPair() :
    p(-1), // point
    p_is_obs(0), // 0 or 1
    q(-1), // face
    q_is_obs(0) // 0 or 1
	{}

void EmbeddedMeshCollision::set_obstacles(
	const float *v0,
	const float *v1,
	int nv,
	const unsigned int *faces,
	int nf)
{
	if (obsdata.V0.rows() != nv)
		obsdata.V0.resize(nv,3);

	if (obsdata.V1.rows() != nv)
		obsdata.V1.resize(nv,3);

	for (int i=0; i<nv; ++i)
	{
		for (int j=0; j<3; ++j)
		{
			obsdata.V0(i,j) = v0[i*3+j];
			obsdata.V1(i,j) = v1[i*3+j];
		}
	}

	if (obsdata.F.rows() != nf)
	{
		obsdata.F.resize(nf,3);
		obsdata.aabbs.resize(nf);
	}

	for (int i=0; i<nf; ++i)
	{
		obsdata.aabbs[i].setEmpty();
		for (int j=0; j<3; ++j)
		{
			int fj = faces[i*3+j];
			obsdata.F(i,j) = fj;
			obsdata.aabbs[i].extend(obsdata.V0.row(fj).transpose());
			obsdata.aabbs[i].extend(obsdata.V1.row(fj).transpose());
		}
	}

	obsdata.tree.init(obsdata.aabbs);

} // end add obstacle

typedef struct DetectThreadData {
	const EmbeddedMeshData *mesh;
	const EmbeddedMeshCollision::ObstacleData *obsdata;
	const Eigen::MatrixXd *x0;
	const Eigen::MatrixXd *x1;
	std::vector<std::vector<VFCollisionPair> > *pt_vf_pairs; // per thread pairs
} DetectThreadData;

static void parallel_detect(
	void *__restrict userdata,
	const int i,
	const TaskParallelTLS *__restrict tls)
{
	// Comments say "don't use this" but how else am I supposed
	// to get the thread ID?
	int thread_idx = BLI_task_parallel_thread_id(tls);
	DetectThreadData *td = (DetectThreadData*)userdata;
	std::vector<VFCollisionPair> &tl_pairs = td->pt_vf_pairs->at(thread_idx);

	int tet_idx = td->mesh->vtx_to_tet[i];
	RowVector4i tet = td->mesh->tets.row(tet_idx);
	Vector4d bary = td->mesh->barys.row(i);
	
	// First, get the surface vertex
	Vector3d pt = 
		bary[0] * td->x1->row(tet[0]) +
		bary[1] * td->x1->row(tet[1]) +
		bary[2] * td->x1->row(tet[2]) +
		bary[3] * td->x1->row(tet[3]);

	// TODO
	// This won't work for overlapping obstacles.
	// We would instead need something like a signed distance field
	// or continuous collision detection.

	PointInTriangleMeshTraverse<double> pt_in_mesh(
		pt, &td->obsdata->V1, &td->obsdata->F);
	td->obsdata->tree.traverse(pt_in_mesh);
	if (pt_in_mesh.output.num_hits() % 2 != 1)
		return;

	// If we are inside an obstacle, we
	// have to project to the nearest surface

} // end parallel lin solve

int EmbeddedMeshCollision::detect(
	const Eigen::MatrixXd *x0,
	const Eigen::MatrixXd *x1)
{
	if (mesh==NULL)
		return 0;

	update_bvh(x0,x1);

	int max_threads = std::max(1,BLI_system_thread_count());
	std::vector<std::vector<VFCollisionPair> > pt_vf_pairs
		(max_threads, std::vector<VFCollisionPair>());

	DetectThreadData thread_data = {
		.mesh = mesh,
		.obsdata = &obsdata,
		.x0 = x0,
		.x1 = x1,
		.pt_vf_pairs = &pt_vf_pairs
	};

	int nv = x1->rows();
	TaskParallelSettings settings;
	BLI_parallel_range_settings_defaults(&settings);
	BLI_task_parallel_range(0, nv, &thread_data, parallel_detect, &settings);

	// Combine thread-local results
	vf_pairs.clear();
	for (int i=0; i<max_threads; ++i)
	{
		const std::vector<VFCollisionPair> &tl_pairs = pt_vf_pairs[i];
		vf_pairs.insert(vf_pairs.end(), tl_pairs.begin(), tl_pairs.end());
	}

	return vf_pairs.size();
}
	/*
void EmbeddedMeshCollision::detect(const Eigen::MatrixXd *x0, const Eigen::MatrixXd *x1)
{

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
		PointInTriangleMeshTraverse<double> pt_in_mesh(pt, &V1, &F);
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

} // end emb collision detect
*/
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
