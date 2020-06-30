// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_collision.h"
#include "admmpd_bvh_traverse.h"
#include "admmpd_math.h"

#include "BLI_assert.h"
#include "BLI_task.h"
#include "BLI_threads.h"

#include <iostream>
#include <sstream>

namespace admmpd {
using namespace Eigen;

VFCollisionPair::VFCollisionPair() :
    p_idx(-1), // point
    p_is_obs(0), // 0 or 1
    q_idx(-1), // face
    q_is_obs(0), // 0 or 1
	pt_on_q(0,0,0)
	{}

void Collision::set_obstacles(
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

void Collision::detect_discrete_vf(
	const Eigen::Vector3d &pt,
	int pt_idx,
	bool pt_is_obs,
	const AABBTree<double,3> *mesh_tree,
	const Eigen::MatrixXd *mesh_x,
	const Eigen::MatrixXi *mesh_tris,
	bool mesh_is_obs,
	double floor_z,
	std::vector<VFCollisionPair> *pairs)
{
	// Special case, check if we are below the floor
	if (pt[2] < floor_z)
	{
		pairs->emplace_back();
		VFCollisionPair &pair = pairs->back();
		pair.p_idx = pt_idx;
		pair.p_is_obs = pt_is_obs;
		pair.q_idx = -1;
		pair.q_is_obs = 1;
		pair.pt_on_q = Vector3d(pt[0],pt[1],floor_z);
	}

	// Any faces to detect against?
	if (mesh_tris->rows()==0)
		return;

	// TODO
	// This won't work for overlapping obstacles.
	// We would instead need something like a signed distance field
	// or continuous collision detection.

	PointInTriangleMeshTraverse<double> pt_in_mesh(
		pt, mesh_x, mesh_tris);
	mesh_tree->traverse(pt_in_mesh);
	if (pt_in_mesh.output.num_hits() % 2 != 1)
		return;

	// If we are inside an obstacle, we
	// have to project to the nearest surface.

	// TODO
	// Consider replacing this with BLI codes:
	// BLI_bvhtree_find_nearest in BLI_kdopbvh.h
	// But since it doesn't have a point-in-mesh
	// detection, we may as use our own BVH/traverser
	// for now.

	NearestTriangleTraverse<double> find_nearest_tri(
		pt, mesh_x, mesh_tris);
	mesh_tree->traverse(find_nearest_tri);
	if (find_nearest_tri.output.prim < 0) // error
		return;

	// Create collision pair
	pairs->emplace_back();
	VFCollisionPair &pair = pairs->back();
	pair.p_idx = pt_idx;
	pair.p_is_obs = pt_is_obs;
	pair.q_idx = find_nearest_tri.output.prim;
	pair.q_is_obs = mesh_is_obs;
	pair.pt_on_q = find_nearest_tri.output.pt_on_tri;

	// Compute face normal
//	BLI_assert(pair.q_idx >= 0);
//	BLI_assert(pair.q_idx < mesh_tris->rows());
//	RowVector3i tri_inds = mesh_tris->row(pair.q_idx);
//	BLI_assert(tri_inds[0]>=0 && tri_inds[0]<mesh_x->rows());
//	BLI_assert(tri_inds[1]>=0 && tri_inds[1]<mesh_x->rows());
//	BLI_assert(tri_inds[2]>=0 && tri_inds[2]<mesh_x->rows());
//	Vector3d tri[3] = {
//		mesh_x->row(tri_inds[0]),
//		mesh_x->row(tri_inds[1]),
//		mesh_x->row(tri_inds[2])
//	};

//	std::stringstream ss;
//	ss << "\nhit:" <<
//		"\n\t normal: " << pair.q_n.transpose() <<
//		"\n\t pt: " << pt.transpose() <<
//		"\n\t pt_on_tri: " << pair.pt_on_q.transpose() <<
//		std::endl;
//	printf("%s",ss.str().c_str());

} // end detect_discrete_vf

typedef struct DetectThreadData {
	const TetMeshData *tetmesh;
	const EmbeddedMeshData *embmesh;
	const Collision::ObstacleData *obsdata;
	const Eigen::MatrixXd *x0;
	const Eigen::MatrixXd *x1;
	double floor_z;
	std::vector<std::vector<VFCollisionPair> > *pt_vf_pairs; // per thread pairs
} DetectThreadData;

static void parallel_detect(
	void *__restrict userdata,
	const int i,
	const TaskParallelTLS *__restrict tls)
{
	DetectThreadData *td = (DetectThreadData*)userdata;

	// Comments say "don't use this" but how else am I supposed
	// to get the thread ID? It appears to return the right thing...
	int thread_idx = BLI_task_parallel_thread_id(tls);
	std::vector<VFCollisionPair> &tl_pairs = td->pt_vf_pairs->at(thread_idx);

	if (td->tetmesh != nullptr)
	{

	} // end detect with tet meshes
	
	if (td->embmesh != nullptr)
	{
		int tet_idx = td->embmesh->vtx_to_tet[i];
		RowVector4i tet = td->embmesh->tets.row(tet_idx);
		Vector4d bary = td->embmesh->barys.row(i);
		Vector3d pt = 
			bary[0] * td->x1->row(tet[0]) +
			bary[1] * td->x1->row(tet[1]) +
			bary[2] * td->x1->row(tet[2]) +
			bary[3] * td->x1->row(tet[3]);

		Collision::detect_discrete_vf(
			pt, i, false,
			&td->obsdata->tree,
			&td->obsdata->V1,
			&td->obsdata->F,
			true,
			td->floor_z,
			&tl_pairs );

	} // end detect with embedded meshes

} // end parallel detect

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
		.tetmesh = nullptr,
		.embmesh = mesh,
		.obsdata = &obsdata,
		.x0 = x0,
		.x1 = x1,
		.floor_z = floor_z,
		.pt_vf_pairs = &pt_vf_pairs
	};

	int nev = mesh->x_rest.rows();

	TaskParallelSettings settings;
	BLI_parallel_range_settings_defaults(&settings);
	BLI_task_parallel_range(0, nev, &thread_data, parallel_detect, &settings);

	// Combine thread-local results
	vf_pairs.clear();
	for (int i=0; i<max_threads; ++i)
	{
		const std::vector<VFCollisionPair> &tl_pairs = pt_vf_pairs[i];
		vf_pairs.insert(vf_pairs.end(), tl_pairs.begin(), tl_pairs.end());
	}

	return vf_pairs.size();
} // end detect

void EmbeddedMeshCollision::linearize(
	const Eigen::MatrixXd *x,
	std::vector<Eigen::Triplet<double> > *trips,
	std::vector<double> *d)
{
	BLI_assert(x != NULL);
	BLI_assert(x->cols() == 3);

	int np = vf_pairs.size();
	if (np==0)
		return;

	d->reserve((int)d->size() + np);
	trips->reserve((int)trips->size() + np*3*4);

	for (int i=0; i<np; ++i)
	{
		VFCollisionPair &pair = vf_pairs[i];
		int emb_p_idx = pair.p_idx;


    //p_idx(-1), // point
    //p_is_obs(0), // 0 or 1
    //q_idx(-1), // face
    //q_is_obs(0), // 0 or 1
//	pt_on_q(0,0,0)

		// Compute normal of triangle at x
		Vector3d q_n(0,0,0);
		RowVector3i q_inds(-1,-1,-1);
		std::vector<Vector3d> q_tris;
		if (pair.q_is_obs)
		{
			// Special case, floor
			if (pair.q_idx == -1)
			{
				q_n = Vector3d(0,0,1);
			}
			else
			{
				q_inds = obsdata.F.row(pair.q_idx);
				q_tris = {
					obsdata.V1.row(q_inds[0]),
					obsdata.V1.row(q_inds[1]),
					obsdata.V1.row(q_inds[2])
				};
				q_n = ((q_tris[1]-q_tris[0]).cross(q_tris[2]-q_tris[0]));
				q_n.normalize();
			}
		}

		// TODO: obstacle point intersecting mesh
		if (pair.p_is_obs)
		{
			continue;
		}

		// Obstacle collision
		if (pair.q_is_obs)
		{
			RowVector4d bary = mesh->barys.row(emb_p_idx);
			int tet_idx = mesh->vtx_to_tet[emb_p_idx];
			RowVector4i tet = mesh->tets.row(tet_idx);
			int c_idx = d->size();
			d->emplace_back(q_n.dot(pair.pt_on_q));
			for (int j=0; j<4; ++j)
			{
				trips->emplace_back(c_idx, tet[j]*3+0, bary[j]*q_n[0]);
				trips->emplace_back(c_idx, tet[j]*3+1, bary[j]*q_n[1]);
				trips->emplace_back(c_idx, tet[j]*3+2, bary[j]*q_n[2]);
			}
			continue;
		}

		// Self collisions
	
	} // end loop pairs

} // end jacobian

/*
int TetMeshCollision::detect(
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
		.tetmesh = mesh,
		.embmesh = nullptr,
		.obsdata = &obsdata,
		.x0 = x0,
		.x1 = x1,
		.floor_z = floor_z,
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

void TetMeshCollision::jacobian(
	const Eigen::MatrixXd *x,
	std::vector<Eigen::Triplet<double> > *trips_x,
	std::vector<Eigen::Triplet<double> > *trips_y,
	std::vector<Eigen::Triplet<double> > *trips_z,
	std::vector<double> *l)
{
	BLI_assert(x != NULL);
	BLI_assert(x->cols() == 3);

	int np = vf_pairs.size();
	if (np==0)
		return;

	l->reserve((int)l->size() + np);
	trips_x->reserve((int)trips_x->size() + np*4);
	trips_y->reserve((int)trips_y->size() + np*4);
	trips_z->reserve((int)trips_z->size() + np*4);

	for (int i=0; i<np; ++i)
	{
		VFCollisionPair &pair = vf_pairs[i];

		// TODO: obstacle point intersecting mesh
		if (pair.p_is_obs)
		{
			continue;
		}

		// Obstacle collision
		if (pair.q_is_obs)
		{
			int c_idx = l->size();
			bool has_qnx = std::abs(pair.q_n[0])>0;
			bool has_qny = std::abs(pair.q_n[1])>0;
			bool has_qnz = std::abs(pair.q_n[2])>0;
			if (has_qnx)
				trips_x->emplace_back(c_idx, pair.p_idx, pair.q_n[0]);
			if (has_qny)
				trips_y->emplace_back(c_idx, pair.p_idx, pair.q_n[1]);
			if (has_qnz)
				trips_z->emplace_back(c_idx, pair.p_idx, pair.q_n[2]);
			l->emplace_back( pair.q_n.dot(pair.pt_on_q));
			continue;
		}
	
	} // end loop pairs

} // end jacobian of tet mesh intersect
*/

} // namespace admmpd
