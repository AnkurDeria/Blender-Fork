


#ifndef _ADMMPD_MATH_H
#define _ADMMPD_MATH_H

#include <Eigen/Geometry>

namespace admmpd {
namespace barycoords {

Eigen::Matrix<double,4,1> point_tet(
    const Eigen::Matrix<double,3,1> &p,
    const Eigen::Matrix<double,3,1> &a,
    const Eigen::Matrix<double,3,1> &b,
    const Eigen::Matrix<double,3,1> &c,
    const Eigen::Matrix<double,3,1> &d);

bool point_in_tet(
    const Eigen::Matrix<double,3,1> &p,
    const Eigen::Matrix<double,3,1> &a,
    const Eigen::Matrix<double,3,1> &b,
    const Eigen::Matrix<double,3,1> &c,
    const Eigen::Matrix<double,3,1> &d);

} // namespace barycoords

} // namespace admmpd

#endif //_ADMMPD_MATH_H
