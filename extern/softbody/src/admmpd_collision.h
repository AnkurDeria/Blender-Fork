// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_COLLISION_H_
#define ADMMPD_COLLISION_H_

#include <Eigen/Sparse>
#include <vector>

namespace admmpd {

// I'll update this class/structure another day.
// For now let's get something in place to do floor collisions.
// Probably will work better to use uber-collision class for
// all self and obstacle collisions, reducing the amount of
// for-all vertices loops.
class Collider {
public:
    virtual void detect(
        const Eigen::MatrixXd *x) = 0;
    virtual void jacobian(
        const Eigen::MatrixXd *x,
    	std::vector<Eigen::Triplet<double> > *trips_x,
        std::vector<Eigen::Triplet<double> > *trips_y,
    	std::vector<Eigen::Triplet<double> > *trips_z,
		std::vector<double> *l) = 0;
};

class FloorCollider : public Collider {
public:
    void detect(const Eigen::MatrixXd *x);
    void jacobian(
        const Eigen::MatrixXd *x,
    	std::vector<Eigen::Triplet<double> > *trips_x,
        std::vector<Eigen::Triplet<double> > *trips_y,
    	std::vector<Eigen::Triplet<double> > *trips_z,
		std::vector<double> *l);
};

} // namespace admmpd

#endif // ADMMPD_COLLISION_H_
