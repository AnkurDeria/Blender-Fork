// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_lattice.h"
#include "admmpd_math.h"
#include "admmpd_bvh.h"
#include "admmpd_bvh_traverse.h"
#include <iostream>
#include <unordered_map>
#include "BLI_task.h" // threading

namespace admmpd {
using namespace Eigen;


bool Lattice::generate(
	const Eigen::MatrixXd &V,
    Eigen::MatrixXd *X, // lattice vertices, n x 3
    Eigen::MatrixXi *T) // lattice elements, m x 4
{
	AlignedBox<double,3> aabb;
	int nv = V.rows();
	for(int i=0; i<nv; ++i)
		aabb.extend(V.row(i).transpose());
	
	aabb.extend(aabb.min()-Vector3d::Ones()*1e-6);
	aabb.extend(aabb.max()+Vector3d::Ones()*1e-6);
	Vector3d mesh_boxmin = aabb.min();
	Vector3d sizes = aabb.sizes();
	double grid_dx = sizes.maxCoeff() * 0.2;

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

	// Copy data into matrices
	{
		nv = verts.size();
		X->resize(nv,3);
		for(int i=0; i<nv; ++i){
			for(int j=0; j<3; ++j){
				X->operator()(i,j) = verts[i][j];
			}
		}
		int nt = tets.size();
		T->resize(nt,4);
		for(int i=0; i<nt; ++i){
			for(int j=0; j<4; ++j){
				T->operator()(i,j) = tets[i][j];
			}
		}
	}

	return compute_vtx_tet_mapping(
		&V, &vtx_to_tet, &barys,
		X, T);

} // end gen lattice

typedef struct FindTetThreadData {
	AABBTree<double,3> *tree;
	const MatrixXd *pts;
	VectorXi *pts_to_tet;
	MatrixXd *barys;
	const MatrixXd *tet_x;
	const MatrixXi *tets;
} FindTetThreadData;

static void parallel_point_in_tet(
	void *__restrict userdata,
	const int i,
	const TaskParallelTLS *__restrict UNUSED(tls))
{
	FindTetThreadData *td = (FindTetThreadData*)userdata;
	Vector3d pt = td->pts->row(i);
	PointInTetTraverse<double> traverser(pt, td->tet_x, td->tets);
	bool success = td->tree->traverse(traverser);
	int tet_idx = traverser.output.prim;
	if (success && tet_idx >= 0)
	{
		RowVector4i tet = td->tets->row(tet_idx);
		Vector3d t[4] = {
			td->tet_x->row(tet[0]),
			td->tet_x->row(tet[1]),
			td->tet_x->row(tet[2]),
			td->tet_x->row(tet[3])
		};
		td->pts_to_tet->operator[](i) = tet_idx;
		Vector4d b = barycoords::point_tet(pt,t[0],t[1],t[2],t[3]);
		td->barys->row(i) = b;
	}
} // end parallel lin solve

bool Lattice::compute_vtx_tet_mapping(
	const Eigen::MatrixXd *vtx_, // embedded vertices, p x 3
	Eigen::VectorXi *vtx_to_tet_, // what tet vtx is embedded in, p x 1
	Eigen::MatrixXd *barys_, // barycoords of the embedding, p x 4
	const Eigen::MatrixXd *x_, // lattice vertices, n x 3
	const Eigen::MatrixXi *tets_) // lattice elements, m x 4
{
	if (!vtx_ || !vtx_to_tet_ || !barys_ || !x_ || !tets_)
		return false;

	int nv = vtx_->rows();
	if (nv==0)
		return false;

	barys_->resize(nv,4);
	barys_->setOnes();
	vtx_to_tet_->resize(nv);
	int nt = tets_->rows();

	// BVH tree for finding point-in-tet and computing
	// barycoords for each embedded vertex
	std::vector<AlignedBox<double,3> > tet_aabbs;
	tet_aabbs.resize(nt);
	for (int i=0; i<nt; ++i)
	{
		tet_aabbs[i].setEmpty();
		RowVector4i tet = tets_->row(i);
		for (int j=0; j<4; ++j)
		{
			tet_aabbs[i].extend(x_->row(tet[j]).transpose());
		}
	}

	AABBTree<double,3> tree;
	tree.init(tet_aabbs);

	FindTetThreadData thread_data = {
		.tree = &tree,
		.pts = vtx_,
		.pts_to_tet = vtx_to_tet_,
		.barys = barys_,
		.tet_x = x_,
		.tets = tets_
	};
	TaskParallelSettings settings;
	BLI_parallel_range_settings_defaults(&settings);
	BLI_task_parallel_range(0, nv, &thread_data, parallel_point_in_tet, &settings);

	// Double check we set (valid) barycoords for every embedded vertex
	const double eps = 1e-8;
	for (int i=0; i<nv; ++i)
	{
		RowVector4d b = barys_->row(i);
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

Eigen::Vector3d Lattice::get_mapped_vertex(
	int idx,
	const Eigen::MatrixXd *x_or_v,
	const Eigen::MatrixXi *tets )
{
    int t_idx = vtx_to_tet[idx];
    RowVector4i tet = tets->row(t_idx);
    RowVector4d b = barys.row(idx);
    return Vector3d(
		x_or_v->row(tet[0]) * b[0] +
		x_or_v->row(tet[1]) * b[1] +
		x_or_v->row(tet[2]) * b[2] +
		x_or_v->row(tet[3]) * b[3]);
}

} // namespace admmpd
