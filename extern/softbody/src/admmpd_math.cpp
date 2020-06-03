

#include "admmpd_math.h"

namespace admmpd {
namespace barycoords {

    Eigen::Matrix<double,4,1> point_tet(
        const Eigen::Vector3d &p,
        const Eigen::Vector3d &a,
        const Eigen::Vector3d &b,
        const Eigen::Vector3d &c,
        const Eigen::Vector3d &d)
    {
		auto scalar_triple_product = [](
			const Eigen::Vector3d &u,
			const Eigen::Vector3d &v,
			const Eigen::Vector3d &w )
		{
			return u.dot(v.cross(w));
		};
		Eigen::Vector3d ap = p - a;
		Eigen::Vector3d bp = p - b;
		Eigen::Vector3d ab = b - a;
		Eigen::Vector3d ac = c - a;
		Eigen::Vector3d ad = d - a;
		Eigen::Vector3d bc = c - b;
		Eigen::Vector3d bd = d - b;
		double va6 = scalar_triple_product(bp, bd, bc);
		double vb6 = scalar_triple_product(ap, ac, ad);
		double vc6 = scalar_triple_product(ap, ad, ab);
		double vd6 = scalar_triple_product(ap, ab, ac);
		double v6 = 1.0 / scalar_triple_product(ab, ac, ad);
		return Eigen::Matrix<double,4,1>(va6*v6, vb6*v6, vc6*v6, vd6*v6);
    } // end point tet barycoords

	// Checks that it's on the "correct" side of the normal
	// for each face of the tet. Assumes winding points inward.
	bool point_in_tet(
		const Eigen::Vector3d &p,
		const Eigen::Vector3d &a,
		const Eigen::Vector3d &b,
		const Eigen::Vector3d &c,
		const Eigen::Vector3d &d)
	{
		using namespace Eigen;
		auto check_face = [](
			const Vector3d &point,
			const Vector3d &p0,
			const Vector3d &p1,
			const Vector3d &p2,
			const Vector3d &p3 )
		{
			Vector3d n = (p1-p0).cross(p2-p0);
			double dp3 = n.dot(p3-p0);
			double dp = n.dot(point-p0);
			return (dp3*dp >= 0);
		};
		return
			check_face(p, a, b, c, d) &&
			check_face(p, b, c, d, a) &&
			check_face(p, c, d, a, b) &&
			check_face(p, d, a, b, c);
	}

} // namespace barycoords
} // namespace admmpd