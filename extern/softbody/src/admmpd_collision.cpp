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
	q_bary(0,0,0)
	{}

void Collision::set_obstacles(
	const float *v0,
	const float *v1,
	int nv,
	const unsigned int *faces,
	int nf)
{
	if (nv==0 || nf==0)
	{
//		obsdata.mesh = admmpd::EmbeddedMesh();
		return;
	}

	if (obsdata.V.rows() != nv)
		obsdata.V.resize(nv,3);
	for (int i=0; i<nv; ++i)
	{
		for (int j=0; j<3; ++j)
			obsdata.V(i,j) = v1[i*3+j];
	}

	if ((int)obsdata.leaves.size() != nf)
		obsdata.leaves.resize(nf);
	if (obsdata.F.rows() != nf)
		obsdata.F.resize(nf,3);
	for (int i=0; i<nf; ++i)
	{
		obsdata.leaves[i].setEmpty();
		for (int j=0; j<3; ++j)
		{
			obsdata.F(i,j) = faces[i*3+j];
			Vector3d vi = obsdata.V.row(obsdata.F(i,j)).transpose();
			obsdata.leaves[i].extend(vi);
		}
	}

	obsdata.tree.init(obsdata.leaves);
/*
  	if (!obsdata.mesh.generate(V,F,true,2))
		return;

	int nt = obsdata.mesh.lat_tets.rows();
	if ((int)obsdata.tet_leaves.size() != nt)
		obsdata.tet_leaves.resize(nt);

	for (int i=0; i<nt; ++i)
	{
		AlignedBox<double,3> &box = obsdata.tet_leaves[i];
		box.setEmpty();
		RowVector4i t = obsdata.mesh.lat_tets.row(i);
		box.extend(obsdata.mesh.lat_rest_x.row(t[0]).transpose());
		box.extend(obsdata.mesh.lat_rest_x.row(t[1]).transpose());
		box.extend(obsdata.mesh.lat_rest_x.row(t[2]).transpose());
		box.extend(obsdata.mesh.lat_rest_x.row(t[3]).transpose());
//		box.extend(box.min()-Vector3d::Ones()*settings.collision_eps);
//		box.extend(box.max()+Vector3d::Ones()*settings.collision_eps);
	}

	obsdata.tet_tree.init(obsdata.tet_leaves);
*/
} // end add obstacle

std::pair<bool,VFCollisionPair>
Collision::detect_against_obs(
        const Eigen::Vector3d &pt,
        const ObstacleData *obs) const
{
	std::pair<bool,VFCollisionPair> ret = 
		std::make_pair(false, VFCollisionPair());

	if (!obs->has_obs())
		return ret;

	PointInTriangleMeshTraverse<double> pt_in_mesh(pt,&obs->V,&obs->F);
	obs->tree.traverse(pt_in_mesh);
	if (pt_in_mesh.output.num_hits()%2==0)
		return ret;

	NearestTriangleTraverse<double> nearest_tri(pt,&obs->V,&obs->F);
	obs->tree.traverse(nearest_tri);

	ret.first = true;
	ret.second.q_idx = nearest_tri.output.prim;
	ret.second.q_is_obs = true;
	ret.second.q_pt = nearest_tri.output.pt_on_tri;
	return ret;
}

int EmbeddedMeshCollision::detect(
	const Eigen::MatrixXd *x0,
	const Eigen::MatrixXd *x1)
{
	if (mesh==NULL)
		return 0;

	// Do we even need to process collisions?
	if (!this->settings.test_floor &&
		!this->settings.self_collision &&
		!obsdata.has_obs())
		return 0;

	// Updates the BVH of deforming mesh
	update_bvh(x0,x1);

	// We store the results of the collisions in a per-vertex buffer.
	// This is a workaround so we can create them in threads.
	int nev = mesh->emb_rest_x.rows();
	if ((int)per_vertex_pairs.size() != nev)
		per_vertex_pairs.resize(nev, std::vector<VFCollisionPair>());

	//
	// Thread data for detection
	//
	typedef struct DetectThreadData {
		const Collision::Settings *settings;
		const Collision *collision;
		const TetMeshData *tetmesh;
		const EmbeddedMesh *embmesh;
		const Collision::ObstacleData *obsdata;
		const Eigen::MatrixXd *x0;
		const Eigen::MatrixXd *x1;
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
		Vector3d pt_t1 = 
			bary[0] * td->x1->row(tet[0]) +
			bary[1] * td->x1->row(tet[1]) +
			bary[2] * td->x1->row(tet[2]) +
			bary[3] * td->x1->row(tet[3]);

		// Special case, check if we are below the floor
		if (td->settings->test_floor)
		{
			if (pt_t1[2] < td->settings->floor_z)
			{
				vi_pairs.emplace_back();
				VFCollisionPair &pair = vi_pairs.back();
				pair.p_idx = vi;
				pair.p_is_obs = false;
				pair.q_idx = -1;
				pair.q_is_obs = 1;
				pair.q_bary.setZero();
				pair.q_pt = Vector3d(pt_t1[0],pt_t1[1],td->settings->floor_z);
			}
		}

		std::pair<bool,VFCollisionPair> pt_hit_obs =
			td->collision->detect_against_obs(pt_t1,td->obsdata);

		if (pt_hit_obs.first)
		{
			pt_hit_obs.second.p_idx = vi;
			pt_hit_obs.second.p_is_obs = false;
			vi_pairs.emplace_back(pt_hit_obs.second);
		}

	}; // end detect for a single embedded vertex

	DetectThreadData thread_data = {
		.settings = &settings,
		.collision = this,
		.tetmesh = nullptr,
		.embmesh = mesh,
		.obsdata = &obsdata,
		.x0 = x0,
		.x1 = x1,
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

		//
		// If we collided with an obstacle
		//
		if (pair.q_is_obs)
		{
			Vector3d q_n(0,0,0);

			// Special case, floor
			if (pair.q_idx == -1)
			{
				q_n = Vector3d(0,0,1);
			}
			else
			{
				RowVector3i q_inds = obsdata.F.row(pair.q_idx);
				Vector3d q_tris[3] = {
					obsdata.V.row(q_inds[0]),
					obsdata.V.row(q_inds[1]),
					obsdata.V.row(q_inds[2])
				};
				q_n = ((q_tris[1]-q_tris[0]).cross(q_tris[2]-q_tris[0]));
				q_n.normalize();
			}

			// Get the four deforming verts that embed
			// the surface vertices, and add constraints on those.
			RowVector4d bary = mesh->emb_barys.row(emb_p_idx);
			int tet_idx = mesh->emb_vtx_to_tet[emb_p_idx];
			RowVector4i tet = mesh->lat_tets.row(tet_idx);
			int c_idx = d->size();
			d->emplace_back(q_n.dot(pair.q_pt));
			for (int j=0; j<4; ++j)
			{
				trips->emplace_back(c_idx, tet[j]*3+0, bary[j]*q_n[0]);
				trips->emplace_back(c_idx, tet[j]*3+1, bary[j]*q_n[1]);
				trips->emplace_back(c_idx, tet[j]*3+2, bary[j]*q_n[2]);
			}
			continue;
		}

	} // end loop pairs

} // end jacobian

} // namespace admmpd

/*
std::pair<bool,VFCollisionPair>
Collision::detect_point_against_mesh(
        int pt_idx,
        bool pt_is_obs,
        const Eigen::Vector3d &pt,
        bool mesh_is_obs,
        const EmbeddedMesh *emb_mesh,
        const Eigen::MatrixXd *mesh_tets_x,
        const AABBTree<double,3> *mesh_tets_tree) const
{
	std::pair<bool,VFCollisionPair> ret = 
		std::make_pair(false, VFCollisionPair());

	if (mesh_tets_x->rows()==0)
		return ret;

	// Point in tet?
	PointInTetMeshTraverse<double> pt_in_tet(pt,mesh_tets_x,&emb_mesh->lat_tets);
	bool in_mesh = mesh_tets_tree->traverse(pt_in_tet);
	if (!in_mesh)
		return ret;

	// Transform to rest shape
	int tet_idx = pt_in_tet.output.prim;
	RowVector4i tet = emb_mesh->lat_tets.row(tet_idx);
	Vector4d barys = geom::point_tet_barys(pt,
		mesh_tets_x->row(tet[0]), mesh_tets_x->row(tet[1]),
		mesh_tets_x->row(tet[2]), mesh_tets_x->row(tet[3]));
	Vector3d rest_pt =
		barys[0]*emb_mesh->lat_rest_x.row(tet[0])+
		barys[1]*emb_mesh->lat_rest_x.row(tet[1])+
		barys[2]*emb_mesh->lat_rest_x.row(tet[2])+
		barys[3]*emb_mesh->lat_rest_x.row(tet[3]);

	// We are inside the lattice. Find nearest
	// face and see if we are on the inside of the mesh.
	NearestTriangleTraverse<double> nearest_tri(rest_pt,
		&emb_mesh->emb_rest_x,&emb_mesh->emb_faces);
	emb_mesh->emb_rest_tree.traverse(nearest_tri);

	if (nearest_tri.output.prim < 0)
		throw std::runtime_error("detect_point_against_mesh failed to project out");

	RowVector3i f = emb_mesh->emb_faces.row(nearest_tri.output.prim);
	Vector3d fv[3] = {
		emb_mesh->emb_rest_x.row(f[0]),
		emb_mesh->emb_rest_x.row(f[1]),
		emb_mesh->emb_rest_x.row(f[2])
	};
	Vector3d n = (fv[1]-fv[0]).cross(fv[2]-fv[0]);
	n.normalize();
	if ((rest_pt-fv[0]).dot(n)>0)
		return ret; // outside of surface

	ret.first = true;
	ret.second.p_idx = pt_idx;
	ret.second.p_is_obs = pt_is_obs;
	ret.second.q_idx = 
	ret.second.q_is_obs = mesh_is_obs;
	ret.second.q_bary = geom::point_triangle_barys<double>(
		nearest_tri.output.pt_on_tri, fv[0], fv[1], fv[2]);

	return ret;
} // end detect_point_against_mesh
*/