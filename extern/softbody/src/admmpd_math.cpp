// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#include "admmpd_math.h"

namespace admmpd {
using namespace Eigen;

namespace barycoords {

    Matrix<double,4,1> point_tet(
        const Vector3d &p,
        const Vector3d &a,
        const Vector3d &b,
        const Vector3d &c,
        const Vector3d &d)
    {
		auto scalar_triple_product = [](
			const Vector3d &u,
			const Vector3d &v,
			const Vector3d &w )
		{
			return u.dot(v.cross(w));
		};
		Vector3d ap = p - a;
		Vector3d bp = p - b;
		Vector3d ab = b - a;
		Vector3d ac = c - a;
		Vector3d ad = d - a;
		Vector3d bc = c - b;
		Vector3d bd = d - b;
		double va6 = scalar_triple_product(bp, bd, bc);
		double vb6 = scalar_triple_product(ap, ac, ad);
		double vc6 = scalar_triple_product(ap, ad, ab);
		double vd6 = scalar_triple_product(ap, ab, ac);
		double v6 = 1.0 / scalar_triple_product(ab, ac, ad);
		return Matrix<double,4,1>(va6*v6, vb6*v6, vc6*v6, vd6*v6);
    } // end point tet barycoords

	// From Real-Time Collision Detection by Christer Ericson
    Vector3d point_triangle(
        const Vector3d &p,
        const Vector3d &a,
        const Vector3d &b,
        const Vector3d &c)
    {
		Vector3d v0 = b - a;
		Vector3d v1 = c - a;
		Vector3d v2 = p - a;
		double d00 = v0.dot(v0);
		double d01 = v0.dot(v1);
		double d11 = v1.dot(v1);
		double d20 = v2.dot(v0);
		double d21 = v2.dot(v1);
		double denom = d00 * d11 - d01 * d01;
		if (std::abs(denom)<=0)
			return Vector3d::Zero();
		Vector3d r;
		r[1] = (d11 * d20 - d01 * d21) / denom;
		r[2] = (d00 * d21 - d01 * d20) / denom;
		r[0] = 1.0 - r[1] - r[2];
		return r;
	} // end point triangle barycoords

} // namespace barycoords

// Checks that it's on the "correct" side of the normal
// for each face of the tet. Assumes winding points inward.
bool point_in_tet(
	const Vector3d &p,
	const Vector3d &a,
	const Vector3d &b,
	const Vector3d &c,
	const Vector3d &d)
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

} // namespace admmpd
