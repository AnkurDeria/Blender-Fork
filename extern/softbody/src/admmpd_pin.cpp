// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_pin.h"
#include "BLI_assert.h"

namespace admmpd {
using namespace Eigen;

void EmbeddedMeshPin::clear()
{
	pin_k.clear();
	pin_pos.clear();
}

void EmbeddedMeshPin::set_pin(
		int idx,
		const Eigen::Vector3d &qi,
		const Eigen::Vector3d &ki_)
{
	if (!mesh)
		return;

	if (idx<0 || idx>=mesh->emb_rest_x.rows())
		return;

	// Clamp
	Vector3d ki = ki_;
	for (int i=0; i<3; ++i)
		ki[i] = std::max(0.0, ki[i]);

	pin_k[idx] = ki;
	pin_pos[idx] = qi;
}

void EmbeddedMeshPin::linearize(
		const Eigen::MatrixXd *x, // not used yet
		std::vector<Eigen::Triplet<double> > *trips,
		std::vector<double> *q)
{

	(void)(x);
	int np = pin_k.size();
	trips->reserve((int)trips->size() + np*3*4);
	q->reserve((int)q->size() + np*3);

	std::unordered_map<int,Eigen::Vector3d>::const_iterator it_k = pin_k.begin();
	for (; it_k != pin_k.end(); ++it_k)
	{
		int emb_idx = it_k->first;
		const Vector3d &qi = pin_pos[emb_idx];
		const Vector3d &ki = it_k->second;

		int tet_idx = mesh->emb_vtx_to_tet[emb_idx];
		RowVector4d bary = mesh->emb_barys.row(emb_idx);
		RowVector4i tet = mesh->lat_tets.row(tet_idx);

		for (int i=0; i<3; ++i)
		{
			int p_idx = q->size();
			q->emplace_back(qi[i]*ki[i]);
			for (int j=0; j<4; ++j)
				trips->emplace_back(p_idx, tet[j]*3+i, bary[j]*ki[i]);
		}
	}

} // end linearize

//bool EmbeddedMeshPin::has_pin(int idx) const
//{
//	return pin_k.count(idx);
//}

} // namespace admmpd
