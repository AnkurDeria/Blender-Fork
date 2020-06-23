// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_COLLISION_H_
#define ADMMPD_COLLISION_H_

#include "admmpd_bvh.h"
#include "admmpd_types.h"

namespace admmpd {

struct VFCollisionPair {
    int p; // point
    int p_is_obs; // 0 or 1
    int q; // face
    int q_is_obs; // 0 or 1
    VFCollisionPair();
    Eigen::Vector3d barys;
};

// I'll update this class/structure another day.
// For now let's get something in place to do floor collisions.
// Probably will work better to use uber-collision class for
// all self and obstacle collisions, reducing the amount of
// for-all vertices loops.
class Collision {
public:
    virtual ~Collision() {}

    // Returns the number of active constraints
    virtual int detect(
        const Eigen::MatrixXd *x0,
        const Eigen::MatrixXd *x1) = 0;

    virtual void set_obstacles(
        const float *v0,
        const float *v1,
        int nv,
        const unsigned int *faces,
        int nf) = 0;

//    virtual void jacobian(
//        const Eigen::MatrixXd *x,
//    	std::vector<Eigen::Triplet<double> > *trips_x,
//        std::vector<Eigen::Triplet<double> > *trips_y,
//    	std::vector<Eigen::Triplet<double> > *trips_z,
//		std::vector<double> *l) = 0;
};

// Collision detection against multiple meshes
class EmbeddedMeshCollision : public Collision {
public:
    EmbeddedMeshCollision(const EmbeddedMeshData *mesh_) :
        mesh(mesh_),
        floor_z(-std::numeric_limits<double>::max())
        {}

    // Obstacle data created in set_obstacles
    struct ObstacleData {
        Eigen::MatrixXd V0, V1;
        Eigen::MatrixXi F;
        std::vector<Eigen::AlignedBox<double,3> > aabbs;
        AABBTree<double,3> tree;
    } obsdata;

    // I don't really like having to switch up interface style, but we'll
    // do so here to avoid copies that would happen in admmpd_api.
    void set_obstacles(
        const float *v0,
        const float *v1,
        int nv,
        const unsigned int *faces,
        int nf);

    // A floor is so common that it makes sense to hard
    // code floor collision instead of using a floor mesh.
    void set_floor(double z)
    {
        floor_z = z;
    };

    // Performs collision detection and stores pairs
    int detect(
        const Eigen::MatrixXd *x0,
        const Eigen::MatrixXd *x1);

    // Linearizes the collision pairs about x
    // for the constraint Kx=l
    void jacobian(
        const Eigen::MatrixXd *x,
    	std::vector<Eigen::Triplet<double> > *trips_x,
        std::vector<Eigen::Triplet<double> > *trips_y,
    	std::vector<Eigen::Triplet<double> > *trips_z,
		std::vector<double> *l)
    {
        
    }

protected:
    // A ptr to the embedded mesh data
    const EmbeddedMeshData *mesh;
    double floor_z;

    // Pairs are compute on detect
    std::vector<VFCollisionPair> vf_pairs;

    // Updates the tetmesh BVH for self collisions.
    // Called by detect()
    // TODO
    void update_bvh(
        const Eigen::MatrixXd *x0,
        const Eigen::MatrixXd *x1)
    { (void)(x0); (void)(x1); }
};

} // namespace admmpd

#endif // ADMMPD_COLLISION_H_
