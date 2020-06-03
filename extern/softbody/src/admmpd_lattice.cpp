


#include "admmpd_lattice.h"
#include "admmpd_math.h"

//#include "vdb.h"

namespace admmpd {
using namespace Eigen;

// We store our deformable data in column major to make
// matrix-vector mults and solves faster, but we want
// to map raw data as row major.
inline void map_to_x3d(const std::vector<Vector3d> &x_vec, Eigen::MatrixXd *x)
{
	int nx = x_vec.size();
	if (nx==0)
	{
		*x = MatrixXd();
		return;
	}

	x->resize(nx,3);
	for (int i=0; i<nx; ++i)
	{
		x->operator()(i,0) = x_vec[i][0];
		x->operator()(i,1) = x_vec[i][1];
		x->operator()(i,2) = x_vec[i][2];
	}
}

inline void map_to_x4i(const std::vector<Vector4i> &x_vec, Eigen::MatrixXi *x)
{
	int nx = x_vec.size();
	if (nx==0)
	{
		*x = MatrixXi();
		return;
	}

	x->resize(nx,4);
	for (int i=0; i<nx; ++i)
	{
		x->operator()(i,0) = x_vec[i][0];
		x->operator()(i,1) = x_vec[i][1];
		x->operator()(i,2) = x_vec[i][2];
		x->operator()(i,3) = x_vec[i][3];
	}
}

bool Lattice::generate(
	const Eigen::MatrixXd &V,
    Eigen::MatrixXd *x, // lattice vertices, n x 3
    Eigen::MatrixXi *tets) // lattice elements, m x 4
{
	// All vertices enclosed in 5 tets!
	// Will generate actual lattice in the future.
	AlignedBox<double,3> box;
	vtx = V;
	int nv = vtx.rows();
	for (int i=0; i<nv; ++i)
	{
		Vector3d v = vtx.row(i);
		box.extend(v);
	}
	box.extend(box.min() - Vector3d::Ones() * 1e-12);
	box.extend(box.max() + Vector3d::Ones() * 1e-12);

	std::vector<Vector3d> x_vec;
	std::vector<Vector4i> t_vec;
	create_packed_tets(box.min(),box.max(),x_vec,t_vec);

	// Copy vector data to output
	map_to_x3d(x_vec, x);
	map_to_x4i(t_vec, tets);
	return compute_vtx_tet_mapping(&vtx, &vtx_to_tet, &barys, x, tets);

} // end gen lattice

//
// Original source (BSD-2):
// github.com/mattoverby/mclscene/blob/master/include/MCL/EmbeddedMesh.hpp
//
void Lattice::create_packed_tets(
	const Eigen::Vector3d &min,
	const Eigen::Vector3d &max,
	std::vector<Vector3d> &verts,
	std::vector<Vector4i> &tets )
{

	// Top plane, clockwise looking down
	Vector3d a = max;
	Vector3d b( min[0], max[1], max[2] );
	Vector3d c( min[0], max[1], min[2] );
	Vector3d d( max[0], max[1], min[2] );

	// Bottom plan, clockwise looking down
	Vector3d e( max[0], min[1], max[2] );
	Vector3d f( min[0], min[1], max[2] );
	Vector3d g( min[0], min[1], min[2] );
	Vector3d h( max[0], min[1], min[2] );

	// Add the verts
	int nv = verts.size();
	verts.emplace_back( a ); // 0
	verts.emplace_back( b ); // 1
	verts.emplace_back( c ); // 2
	verts.emplace_back( d ); // 3
	verts.emplace_back( e ); // 4
	verts.emplace_back( f ); // 5
	verts.emplace_back( g ); // 6
	verts.emplace_back( h ); // 7

	// Pack 5 tets into the cube
	Vector4i t0( 0, 5, 7, 4 );
	Vector4i t1( 5, 7, 2, 0 );
	Vector4i t2( 5, 0, 2, 1 );
	Vector4i t3( 7, 2, 0, 3 );
	Vector4i t4( 5, 2, 7, 6 );
	Vector4i offset(nv,nv,nv,nv);

	// Add the tets
	tets.emplace_back( t0+offset );
	tets.emplace_back( t1+offset );
	tets.emplace_back( t2+offset );
	tets.emplace_back( t3+offset );
	tets.emplace_back( t4+offset );

//vdb_point(a[0],a[1],a[2]);
//vdb_point(b[0],b[1],b[2]);
//vdb_point(c[0],c[1],c[2]);
//vdb_point(d[0],d[1],d[2]);
//vdb_point(e[0],e[1],e[2]);
//vdb_point(f[0],f[1],f[2]);
//vdb_point(g[0],g[1],g[2]);
//vdb_point(h[0],h[1],h[2]);
//vdb_line(float x0, float y0, float z0, 
//             float x1, float y1, float z1);

} // end create packed tets

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

	// TODO:
	// Use AABB Tree to do point-in-tet and compute bary weighting
	// For now the dumb approach and loop all.

	barys_->resize(nv,4);
	barys_->setZero();
	vtx_to_tet_->resize(nv);
	int nt = tets_->rows();

	for (int i=0; i<nv; ++i)
	{
		Vector3d v = vtx_->row(i);
		for (int j=0; j<nt; ++j)
		{
			RowVector4i tet = tets_->row(j);
			Vector3d t[4] = {
				x_->row(tet[0]),
				x_->row(tet[1]),
				x_->row(tet[2]),
				x_->row(tet[3])
			};
			if (!barycoords::point_in_tet(v,t[0],t[1],t[2],t[3]))
				continue;

			Vector4d b = barycoords::point_tet(v,t[0],t[1],t[2],t[3]);
			vtx_to_tet_->operator[](i) = j;
			barys_->row(i) = b;
			break;
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
/*
void Lattice::map_to_object(
	const Eigen::MatrixXd *x,
	const Eigen::MatrixXi *tets,
	float (*vertexCos)[3])
{
  int nv = vtx.rows();
  for (int i=0; i<nv; ++i)
  {
    int t_idx = vtx_to_tet[i];
    RowVector4i tet = tets->row(t_idx);
    RowVector4d b = barys.row(i);
    vtx.row(i) =
		x->row(tet[0]) * b[0] +
		x->row(tet[1]) * b[1] +
		x->row(tet[2]) * b[2] +
		x->row(tet[3]) * b[3];
	vertexCos[i][0] = vtx(i,0);
	vertexCos[i][1] = vtx(i,1);
	vertexCos[i][2] = vtx(i,2);
  }
} // end map to object
*/
} // namespace admmpd