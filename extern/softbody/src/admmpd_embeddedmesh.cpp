// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_embeddedmesh.h"
#include "admmpd_math.h"
#include "admmpd_bvh.h"
#include "admmpd_bvh_traverse.h"
#include <iostream>
#include <unordered_map>
#include <set>
#include "BLI_task.h" // threading
#include "BLI_assert.h"

namespace admmpd {
using namespace Eigen;

typedef struct KeepTetThreadData {
	const AABBTree<double,3> *tree; // of embedded faces
	const MatrixXd *pts; // of embedded verts
	const MatrixXi *faces; // embedded faces
	const std::vector<Vector3d> *tet_x;
	const std::vector<Vector4i> *tets;
	std::vector<int> *keep_tet; // 0 or 1
} KeepTetThreadData;

static void parallel_keep_tet(
	void *__restrict userdata,
	const int i,
	const TaskParallelTLS *__restrict UNUSED(tls))
{
	KeepTetThreadData *td = (KeepTetThreadData*)userdata;
	RowVector4i tet = td->tets->at(i);
	Vector3d tet_pts[4] = {
		td->tet_x->at(tet[0]),
		td->tet_x->at(tet[1]),
		td->tet_x->at(tet[2]),
		td->tet_x->at(tet[3])
	};

	// Returns true if the tet intersects the
	// surface mesh. Even if it doesn't, we want to keep
	// the tet if it's totally enclosed in the mesh.
	TetIntersectsMeshTraverse<double> tet_hits_mesh(
		tet_pts, td->pts, td->faces);
	bool hit = td->tree->traverse(tet_hits_mesh);
	if (!hit)
	{
		// We only need to check if one vertex of the
		// tet is inside the mesh. If a subset of
		// vertices are inside the mesh, then there would
		// be tri/tri intersection.
		PointInTriangleMeshTraverse<double> pt_in_mesh(
			tet_pts[0], td->pts, td->faces);
		td->tree->traverse(pt_in_mesh);
		if (pt_in_mesh.output.num_hits() % 2 == 1)
		{
			hit = true;
		}
	}

	if (hit) { td->keep_tet->at(i) = 1; }
	else { td->keep_tet->at(i) = 0; }

} // end parallel test if keep tet

bool EmbeddedMesh::generate(
	const Eigen::MatrixXd &V, // embedded verts
	const Eigen::MatrixXi &F, // embedded faces
	EmbeddedMeshData *emb_mesh, // where embedding is stored
	Eigen::MatrixXd *x_tets) // lattice vertices, n x 3
{
	// How big the grid cells are as a fraction
	// of the total mesh.
	static const double GRID_FRAC = 0.2;

	if (emb_mesh==NULL)
		return false;
	if (x_tets==NULL)
		return false;

	emb_mesh->faces = F;
	emb_mesh->x_rest = V;

	AlignedBox<double,3> aabb;
	int nv = V.rows();
	for(int i=0; i<nv; ++i)
		aabb.extend(V.row(i).transpose());
	
	aabb.extend(aabb.min()-Vector3d::Ones()*1e-6);
	aabb.extend(aabb.max()+Vector3d::Ones()*1e-6);
	Vector3d mesh_boxmin = aabb.min();
	Vector3d sizes = aabb.sizes();
	double grid_dx = sizes.maxCoeff() * GRID_FRAC;

	// Generate vertices and tets
	std::vector<Vector3d> verts;
	std::vector<Vector4i> tets;
	{
		std::unordered_map<std::string, int> vertex_map; // [i,j,k]->index

		auto grid_hash = [&](const Vector3d &x)
		{
			Vector3i ll = Vector3i( // nearest gridcell
				std::round((x[0]-mesh_boxmin[0])/grid_dx),
				std::round((x[1]-mesh_boxmin[1])/grid_dx),
				std::round((x[2]-mesh_boxmin[2])/grid_dx));
			return std::to_string(ll[0])+' '+std::to_string(ll[1])+' '+std::to_string(ll[2]);
		};

		auto add_box = [&](int i_, int j_, int k_)
		{
			Vector3d min = mesh_boxmin + Vector3d(i_*grid_dx, j_*grid_dx, k_*grid_dx);
			Vector3d max = mesh_boxmin + Vector3d((i_+1)*grid_dx, (j_+1)*grid_dx, (k_+1)*grid_dx);
			std::vector<Vector3d> v = {
				// Top plane, clockwise looking down
				max,
				Vector3d(min[0], max[1], max[2]),
				Vector3d(min[0], max[1], min[2]),
				Vector3d(max[0], max[1], min[2]),
				// Bottom plan, clockwise looking down
				Vector3d(max[0], min[1], max[2]),
				Vector3d(min[0], min[1], max[2]),
				min,
				Vector3d(max[0], min[1], min[2])
			};
			// Add vertices and get indices of the box
			std::vector<int> b;
			for(int i=0; i<8; ++i)
			{
				std::string hash = grid_hash(v[i]);
				if( vertex_map.count(hash)==0 )
				{
					vertex_map[hash] = verts.size();
					verts.emplace_back(v[i]);
				}
				b.emplace_back(vertex_map[hash]);
			}
			// From the box, create five new tets
			std::vector<Vector4i> new_tets = {
				Vector4i( b[0], b[5], b[7], b[4] ),
				Vector4i( b[5], b[7], b[2], b[0] ),
				Vector4i( b[5], b[0], b[2], b[1] ),
				Vector4i( b[7], b[2], b[0], b[3] ),
				Vector4i( b[5], b[2], b[7], b[6] )
			};
			for(int i=0; i<5; ++i)
				tets.emplace_back(new_tets[i]);
		};

		int ni = std::ceil(sizes[0]/grid_dx);
		int nj = std::ceil(sizes[1]/grid_dx);
		int nk = std::ceil(sizes[2]/grid_dx);
		for(int i=0; i<ni; ++i)
		{
			for(int j=0; j<nj; ++j)
			{
				for(int k=0; k<nk; ++k)
				{
					add_box(i,j,k);
				}
			}
		}

	} // end create boxes and vertices

	// We only want to keep tets that are either
	// a) intersecting the surface mesh
	// b) totally inside the surface mesh
	std::set<int> refd_verts;
	{
		int nf = F.rows();
		std::vector<AlignedBox<double,3> > face_aabb(nf);
		for (int i=0; i<nf; ++i)
		{
			RowVector3i f = F.row(i);
			face_aabb[i].setEmpty();
			for (int j=0; j<3; ++j)
				face_aabb[i].extend(V.row(f[j]).transpose());
		}

		int nt0 = tets.size();
		std::vector<int> keep_tet(nt0,1);

		AABBTree<double,3> mesh_tree;
		mesh_tree.init(face_aabb);
		KeepTetThreadData thread_data = {
			.tree = &mesh_tree,
			.pts = &V,
			.faces = &F,
			.tet_x = &verts,
			.tets = &tets,
			.keep_tet = &keep_tet
		};
		TaskParallelSettings settings;
		BLI_parallel_range_settings_defaults(&settings);
		BLI_task_parallel_range(0, nt0, &thread_data, parallel_keep_tet, &settings);

		// Loop over tets and remove as needed.
		// Mark referenced vertices to compute a mapping.
		int tet_idx = 0;
		for (std::vector<Vector4i>::iterator it = tets.begin(); it != tets.end(); ++tet_idx)
		{
			bool keep = keep_tet[tet_idx];
			if (keep)
			{
				const Vector4i &t = *it;
				refd_verts.emplace(t[0]);
				refd_verts.emplace(t[1]);
				refd_verts.emplace(t[2]);
				refd_verts.emplace(t[3]);
				++it;
			}
			else { it = tets.erase(it); }
		}

	} // end remove unnecessary tets

	// Copy data into matrices and remove unreferenced
	{
		// Computing a mapping of vertices from old to new
		// Delete any unreferenced vertices
		std::unordered_map<int,int> vtx_old_to_new;
		int ntv0 = verts.size(); // original num verts
		int ntv1 = refd_verts.size(); // reduced num verts
		BLI_assert(ntv1 <= ntv0);
		x_tets->resize(ntv1,3);
		int vtx_count = 0;
		for (int i=0; i<ntv0; ++i)
		{
			if (refd_verts.count(i)>0)
			{
				for(int j=0; j<3; ++j){
					x_tets->operator()(vtx_count,j) = verts[i][j];
				}
				vtx_old_to_new[i] = vtx_count;
				vtx_count++;
			}
		}

		// Copy tets to matrix data and update vertices
		int nt = tets.size();
		emb_mesh->tets.resize(nt,4);
		for(int i=0; i<nt; ++i){
			for(int j=0; j<4; ++j){
				int old_idx = tets[i][j];
				BLI_assert(vtx_old_to_new.count(old_idx)>0);
				emb_mesh->tets(i,j) = vtx_old_to_new[old_idx];
			}
		}
	}

	// Now compute the baryweighting for embedded vertices
	return compute_embedding(
		emb_mesh, &V, x_tets);

} // end gen lattice

void EmbeddedMesh::compute_masses(
	EmbeddedMeshData *emb_mesh, // where embedding is stored
	const Eigen::MatrixXd *x_embed, // embedded vertices, p x 3
	const Eigen::MatrixXd *x_tets, // lattice vertices, n x 3
	Eigen::VectorXd *masses_tets, // masses of the lattice verts
	double density_kgm3)
{
	BLI_assert(emb_mesh != NULL);
	BLI_assert(x_embed != NULL);
	BLI_assert(x_tets != NULL);
	BLI_assert(x_tets->rows() > 0);
	BLI_assert(x_tets->cols() == 3);
	BLI_assert(masses_tets != NULL);
	BLI_assert(density_kgm3 > 0);

	// TODO
	// map the area of the surface to the tet vertices

	// Source: https://github.com/mattoverby/mclscene/blob/master/include/MCL/TetMesh.hpp
	// Computes volume-weighted masses for each vertex
	// density_kgm3 is the unit-volume density
	int nx = x_tets->rows();
	masses_tets->resize(nx);
	masses_tets->setZero();
	int n_tets = emb_mesh->tets.rows();
	for (int t=0; t<n_tets; ++t)
	{
		RowVector4i tet = emb_mesh->tets.row(t);
		RowVector3d tet_v0 = x_tets->row(tet[0]);
		Matrix3d edges;
		edges.col(0) = x_tets->row(tet[1]) - tet_v0;
		edges.col(1) = x_tets->row(tet[2]) - tet_v0;
		edges.col(2) = x_tets->row(tet[3]) - tet_v0;
		double vol = std::abs((edges).determinant()/6.f);
		double tet_mass = density_kgm3 * vol;
		masses_tets->operator[](tet[0]) += tet_mass / 4.f;
		masses_tets->operator[](tet[1]) += tet_mass / 4.f;
		masses_tets->operator[](tet[2]) += tet_mass / 4.f;
		masses_tets->operator[](tet[3]) += tet_mass / 4.f;
	}

	// Verify masses
	for (int i=0; i<nx; ++i)
	{
		if (masses_tets->operator[](i) <= 0.0)
		{
			printf("**EmbeddedMesh::compute_masses Error: unreferenced vertex\n");
			masses_tets->operator[](i)=1;
		}
	}
} // end compute masses

typedef struct FindTetThreadData {
	AABBTree<double,3> *tree;
	EmbeddedMeshData *emb_mesh; // thread sets vtx_to_tet and barys
	const Eigen::MatrixXd *x_embed;
	const Eigen::MatrixXd *x_tets;
} FindTetThreadData;

static void parallel_point_in_tet(
	void *__restrict userdata,
	const int i,
	const TaskParallelTLS *__restrict UNUSED(tls))
{
	FindTetThreadData *td = (FindTetThreadData*)userdata;
	Vector3d pt = td->x_embed->row(i);
	PointInTetMeshTraverse<double> traverser(pt, td->x_tets, &td->emb_mesh->tets);
	bool success = td->tree->traverse(traverser);
	int tet_idx = traverser.output.prim;
	if (success && tet_idx >= 0)
	{
		RowVector4i tet = td->emb_mesh->tets.row(tet_idx);
		Vector3d t[4] = {
			td->x_tets->row(tet[0]),
			td->x_tets->row(tet[1]),
			td->x_tets->row(tet[2]),
			td->x_tets->row(tet[3])
		};
		td->emb_mesh->vtx_to_tet[i] = tet_idx;
		Vector4d b = barycoords::point_tet(pt,t[0],t[1],t[2],t[3]);
		td->emb_mesh->barys.row(i) = b;
	}
} // end parallel lin solve

bool EmbeddedMesh::compute_embedding(
	EmbeddedMeshData *emb_mesh, // where embedding is stored
	const Eigen::MatrixXd *x_embed, // embedded vertices, p x 3
	const Eigen::MatrixXd *x_tets) // lattice vertices, n x 3
{
	BLI_assert(emb_mesh!=NULL);
	BLI_assert(x_embed!=NULL);
	BLI_assert(x_tets!=NULL);

	int nv = x_embed->rows();
	if (nv==0)
		return false;

	emb_mesh->barys.resize(nv,4);
	emb_mesh->barys.setOnes();
	emb_mesh->vtx_to_tet.resize(nv);
	int nt = emb_mesh->tets.rows();

	// BVH tree for finding point-in-tet and computing
	// barycoords for each embedded vertex
	std::vector<AlignedBox<double,3> > tet_aabbs;
	tet_aabbs.resize(nt);
	Vector3d veta = Vector3d::Ones()*1e-12;
	for (int i=0; i<nt; ++i)
	{
		tet_aabbs[i].setEmpty();
		RowVector4i tet = emb_mesh->tets.row(i);
		for (int j=0; j<4; ++j)
			tet_aabbs[i].extend(x_tets->row(tet[j]).transpose());

		tet_aabbs[i].extend(tet_aabbs[i].min()-veta);
		tet_aabbs[i].extend(tet_aabbs[i].max()+veta);
	}

	AABBTree<double,3> tree;
	tree.init(tet_aabbs);

	FindTetThreadData thread_data = {
		.tree = &tree,
		.emb_mesh = emb_mesh,
		.x_embed = x_embed,
		.x_tets = x_tets
	};
	TaskParallelSettings settings;
	BLI_parallel_range_settings_defaults(&settings);
	BLI_task_parallel_range(0, nv, &thread_data, parallel_point_in_tet, &settings);

	// Double check we set (valid) barycoords for every embedded vertex
	const double eps = 1e-8;
	for (int i=0; i<nv; ++i)
	{
		RowVector4d b = emb_mesh->barys.row(i);
		if (b.minCoeff() < -eps)
		{
			printf("**Lattice::generate Error: negative barycoords\n");
			return false;
		}
		if (b.maxCoeff() > 1 + eps)
		{
			printf("**Lattice::generate Error: max barycoord > 1\n");
			return false;
		}
		if (b.sum() > 1 + eps)
		{
			printf("**Lattice::generate Error: barycoord sum > 1\n");
			return false;
		}
	}

	return true;

} // end compute vtx to tet mapping

Eigen::Vector3d EmbeddedMesh::get_mapped_vertex(
	const EmbeddedMeshData *emb_mesh,
	const Eigen::MatrixXd *x_data,
	int idx)
{
    int t_idx = emb_mesh->vtx_to_tet[idx];
    RowVector4i tet = emb_mesh->tets.row(t_idx);
    RowVector4d b = emb_mesh->barys.row(idx);
    return Vector3d(
		x_data->row(tet[0]) * b[0] +
		x_data->row(tet[1]) * b[1] +
		x_data->row(tet[2]) * b[2] +
		x_data->row(tet[3]) * b[3]);
}

} // namespace admmpd
