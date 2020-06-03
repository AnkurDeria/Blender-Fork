

#include "admmpd_collision.h"

namespace admmpd {

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
	int nx = x->rows();
	for (int i=0; i<nx; ++i)
	{
		Eigen::Vector3d p = x->row(i);
		if (p[2]>0)
			continue;
	
		int c_idx = l->size();
		trips_z->emplace_back(c_idx,i,1.0);
		l->emplace_back(0);
	}
} // end floor collider Jacobian

} // namespace admmpd