


#ifndef _ADMMPD_MATH_H
#define _ADMMPD_MATH_H

#include <Eigen/Geometry>

namespace admmpd {
namespace barycoords {

Eigen::Vector4d point_tet(
    const Eigen::Vector3d &p,
    const Eigen::Vector3d &a,
    const Eigen::Vector3d &b,
    const Eigen::Vector3d &c,
    const Eigen::Vector3d &d);

bool point_in_tet(
    const Eigen::Vector3d &p,
    const Eigen::Vector3d &a,
    const Eigen::Vector3d &b,
    const Eigen::Vector3d &c,
    const Eigen::Vector3d &d);

} // namespace barycoords

} // namespace admmpd

#endif //_ADMMPD_MATH_H
