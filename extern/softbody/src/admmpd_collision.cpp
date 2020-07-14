// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_collision.h"
#include "admmpd_bvh_traverse.h"
#include "admmpd_geom.h"

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

	Vector3d v_eta = Vector3d::Ones()*settings.collision_eps;
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
		obsdata.aabbs[i].extend(obsdata.aabbs[i].min()-v_eta);
		obsdata.aabbs[i].extend(obsdata.aabbs[i].max()+v_eta);
	}

	obsdata.tree.init(obsdata.aabbs);

} // end add obstacle

std::pair<bool,VFCollisionPair>
Collision::detect_point_against_mesh(
        int pt_idx,
        bool pt_is_obs,
        const Eigen::Vector3d &pt_t0,
        const Eigen::Vector3d &pt_t1,
        bool mesh_is_obs,
        const Eigen::MatrixXd *mesh_x,
        const Eigen::MatrixXi *mesh_tris,
        const AABBTree<double,3> *mesh_tree) const
{
	std::pair<bool,VFCollisionPair> ret = 
		std::make_pair(false, VFCollisionPair());

	// Any faces to detect against?
	if (mesh_tris->rows()==0)
		return ret;	

/*
	VFCollisionPair &pair = ret.second;
	pair.p_idx = pt_idx;
	pair.p_is_obs = false;
	pair.q_idx = 0;
	pair.q_is_obs = 1;
	double max_z = mesh_x->col(2).maxCoeff();
	pair.pt_on_q = Vector3d(pt_t1[0],pt_t1[1],max_z);
*/
	return ret;
/*
	// Use a ray that is pos at t0 to t1
	// with t_max being the displacment. If not moving,
	// set the "eps" variable used by traversal which
	// does point-triangle distance, and considers it a hit
	// if the distance is less than the eps.
	Vector3d dx = (pt_t1-pt_t0);
	double t_max = dx.norm();
	double eps = -1; // used as point-triangle distance query
	//if (t_max < settings.collision_eps)
	if (false)
	{
		dx = Vector3d(0,0,1);
		t_max = settings.collision_eps;
		eps = settings.collision_eps;
	}
	dx.normalize();

	// Traverse the BVH
	RayClosestHit<double> ray_hit_mesh(
		pt_t0, dx,
		mesh_x, mesh_tris,
		eps, 0, t_max);
	mesh_tree->traverse(ray_hit_mesh);

//	if (pt_idx==0)
//	{
//		std::cout << "\n\nV0 (z): " << pt_t0[2] << std::endl;
//		std::cout << "dir: " << dx.transpose() << std::endl;
//		std::cout << ray_hit_mesh.output.prim << std::endl;
//		if(ray_hit_mesh.output.prim >= 0)
//			throw std::runtime_error("YAY HIT");
//	}

	if (ray_hit_mesh.output.prim < 0)
		return ret;

	ret.first = true;
	VFCollisionPair &pair = ret.second;
	pair.p_idx = pt_idx;
	pair.p_is_obs = pt_is_obs;
	pair.q_idx = ray_hit_mesh.output.prim;
	pair.q_is_obs = mesh_is_obs;
	RowVector3i tris = mesh_tris->row(pair.q_idx);
	pair.pt_on_q =
		mesh_x->row(tris[0])*ray_hit_mesh.output.barys[0] +
		mesh_x->row(tris[1])*ray_hit_mesh.output.barys[1] +
		mesh_x->row(tris[2])*ray_hit_mesh.output.barys[2];
	return ret;
*/
} // end detect_point_against_mesh

int EmbeddedMeshCollision::detect(
	const Eigen::MatrixXd *x0,
	const Eigen::MatrixXd *x1)
{
	if (mesh==NULL)
		return 0;

	// Do we even need to process collisions?
	if (!this->settings.test_floor &&
		!this->settings.self_collision &&
		this->obsdata.F.rows()==0)
		return 0;

	// Updates the BVH of self-intersection mesh
	update_bvh(x0,x1);

	// We store the results of the collisions in a
	// per-vertex buffer.
	int nev = mesh->emb_rest_x.rows();
	if ((int)per_vertex_pairs.size() != nev)
		per_vertex_pairs.resize(nev, std::vector<VFCollisionPair>());

	// Set the floor z to inf if bool detect says to
	double floor_z = this->settings.floor_z;
	if (!this->settings.test_floor)
		floor_z = -std::numeric_limits<double>::max();

	//
	// Thread data for detection
	//
	typedef struct DetectThreadData {
		const Collision *collision;
		const TetMeshData *tetmesh;
		const EmbeddedMesh *embmesh;
		const Collision::ObstacleData *obsdata;
		const Eigen::MatrixXd *x0;
		const Eigen::MatrixXd *x1;
		double floor_z;
		std::vector<std::vector<VFCollisionPair> > *per_vertex_pairs;
	} DetectThreadData;

	//
	// Detection function for a single embedded vertex
	//
	auto per_embedded_vertex_detect = [](
		void *__restrict userdata,
		const int vi,
		const TaskParallelTLS *__restrict tls)->void
	{
		(void)(tls);
		DetectThreadData *td = (DetectThreadData*)userdata;
		if (td->embmesh == nullptr)
			return;

		std::vector<VFCollisionPair> &vi_pairs = td->per_vertex_pairs->at(vi);
		vi_pairs.clear();

		int tet_idx = td->embmesh->emb_vtx_to_tet[vi];
		RowVector4i tet = td->embmesh->lat_tets.row(tet_idx);
		Vector4d bary = td->embmesh->emb_barys.row(vi);
		Vector3d pt_t0 = 
			bary[0] * td->x0->row(tet[0]) +
			bary[1] * td->x0->row(tet[1]) +
			bary[2] * td->x0->row(tet[2]) +
			bary[3] * td->x0->row(tet[3]);
		Vector3d pt_t1 = 
			bary[0] * td->x1->row(tet[0]) +
			bary[1] * td->x1->row(tet[1]) +
			bary[2] * td->x1->row(tet[2]) +
			bary[3] * td->x1->row(tet[3]);

		// Special case, check if we are below the floor
		if (pt_t1[2] < td->floor_z)
		{
			vi_pairs.emplace_back();
			VFCollisionPair &pair = vi_pairs.back();
			pair.p_idx = vi;
			pair.p_is_obs = false;
			pair.q_idx = -1;
			pair.q_is_obs = 1;
			pair.pt_on_q = Vector3d(pt_t1[0],pt_t1[1],td->floor_z);
		}

		std::pair<bool,VFCollisionPair> pt_hit_obs =
			td->collision->detect_point_against_mesh(vi, false, pt_t0, pt_t1,
			true, &td->obsdata->V1, &td->obsdata->F, &td->obsdata->tree);

		if (pt_hit_obs.first)
			vi_pairs.emplace_back(pt_hit_obs.second);

	}; // end detect for a single embedded vertex

	DetectThreadData thread_data = {
		.collision = this,
		.tetmesh = nullptr,
		.embmesh = mesh,
		.obsdata = &obsdata,
		.x0 = x0,
		.x1 = x1,
		.floor_z = floor_z,
		.per_vertex_pairs = &per_vertex_pairs
	};

	TaskParallelSettings thrd_settings;
	BLI_parallel_range_settings_defaults(&thrd_settings);
	BLI_task_parallel_range(0, nev, &thread_data, per_embedded_vertex_detect, &thrd_settings);

	vf_pairs.clear();
	for (int i=0; i<nev; ++i)
	{
		int pvp = per_vertex_pairs[i].size();
		for (int j=0; j<pvp; ++j)
			vf_pairs.emplace_back(Vector2i(i,j));
	}

	return vf_pairs.size();
} // end detect


void EmbeddedMeshCollision::update_bvh(
	const Eigen::MatrixXd *x0,
	const Eigen::MatrixXd *x1)
{
	(void)(x0);
	(void)(x1);
}


void EmbeddedMeshCollision::graph(
	std::vector<std::set<int> > &g)
{
	int np = vf_pairs.size();
	if (np==0)
		return;

	int nv = mesh->lat_rest_x.rows();
	if ((int)g.size() < nv)
		g.resize(nv, std::set<int>());

	for (int i=0; i<np; ++i)
	{
		Vector2i pair_idx = vf_pairs[i];
		VFCollisionPair &pair = per_vertex_pairs[pair_idx[0]][pair_idx[1]];
		std::set<int> stencil;

		if (!pair.p_is_obs)
		{
			int tet_idx = mesh->emb_vtx_to_tet[pair.p_idx];
			RowVector4i tet = mesh->lat_tets.row(tet_idx);
			stencil.emplace(tet[0]);
			stencil.emplace(tet[1]);
			stencil.emplace(tet[2]);
			stencil.emplace(tet[3]);
		}
		if (!pair.q_is_obs)
		{
			RowVector3i emb_face = mesh->emb_faces.row(pair.q_idx);
			for (int j=0; j<3; ++j)
			{
				int tet_idx = mesh->emb_vtx_to_tet[emb_face[j]];
				RowVector4i tet = mesh->lat_tets.row(tet_idx);
				stencil.emplace(tet[0]);
				stencil.emplace(tet[1]);
				stencil.emplace(tet[2]);
				stencil.emplace(tet[3]);	
			}
		}

		for (std::set<int>::iterator it = stencil.begin();
			it != stencil.end(); ++it)
		{
			for (std::set<int>::iterator it2 = stencil.begin();
				it2 != stencil.end(); ++it2)
			{
				if (*it == *it2)
					continue;
				g[*it].emplace(*it2);
			}
		}
	}
} // end graph

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

	//int nx = x->rows();
	d->reserve((int)d->size() + np);
	trips->reserve((int)trips->size() + np*3*4);

	for (int i=0; i<np; ++i)
	{
		Vector2i pair_idx = vf_pairs[i];
		VFCollisionPair &pair = per_vertex_pairs[pair_idx[0]][pair_idx[1]];
		int emb_p_idx = pair.p_idx;

		// Compute normal of triangle at x
		if (pair.q_is_obs)
		{
			Vector3d q_n(0,0,0);
			RowVector3i q_inds(-1,-1,-1);
			std::vector<Vector3d> q_tris;

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

			RowVector4d bary = mesh->emb_barys.row(emb_p_idx);
			int tet_idx = mesh->emb_vtx_to_tet[emb_p_idx];
			RowVector4i tet = mesh->lat_tets.row(tet_idx);
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

		// TODO: obstacle point intersecting mesh
		if (pair.p_is_obs)
		{
			continue;
		}

		// Self collisions
	
	} // end loop pairs

} // end jacobian

} // namespace admmpd
