// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_COLLISION_H_
#define ADMMPD_COLLISION_H_

#include "admmpd_bvh.h"
#include "admmpd_types.h"

namespace admmpd {

struct VFCollisionPair {
    int p_idx; // point
    int p_is_obs; // 0 or 1
    int q_idx; // face
    int q_is_obs; // 0 or 1
    Eigen::Vector3d pt_on_q; // point of collision on q
//    int floor; // 0 or 1, special case
//    Eigen::Vector3d barys;
    Eigen::Vector3d q_n; // face normal
    VFCollisionPair();
};

class Collision {
public:
    // Obstacle data created in set_obstacles
    struct ObstacleData {
        Eigen::MatrixXd V0, V1;
        Eigen::MatrixXi F;
        std::vector<Eigen::AlignedBox<double,3> > aabbs;
        AABBTree<double,3> tree;
    } obsdata;

    virtual ~Collision() {}

    // Performs collision detection.
    // Returns the number of active constraints.
    virtual int detect(
        const Eigen::MatrixXd *x0,
        const Eigen::MatrixXd *x1) = 0;

    // Set the soup of obstacles for this time step.
    // I don't really like having to switch up interface style, but we'll
    // do so here to avoid copies that would happen in admmpd_api.
    virtual void set_obstacles(
        const float *v0,
        const float *v1,
        int nv,
        const unsigned int *faces,
        int nf);

    // Special case for floor since it's common.
    virtual void set_floor(double z) = 0;

    // Linearize the constraints and return Jacobian tensor.
    // Constraints are linearized about x for constraint
    // K x = l
    virtual void jacobian(
        const Eigen::MatrixXd *x,
    	std::vector<Eigen::Triplet<double> > *trips_x,
        std::vector<Eigen::Triplet<double> > *trips_y,
    	std::vector<Eigen::Triplet<double> > *trips_z,
		std::vector<double> *l) = 0;

    // Given a point and a surface mesh,
    // perform discrete collision and create
    // a vertex-face collision pair if colliding.
    // Also adds collision pairs if below floor.
    static void detect_discrete_vf(
        const Eigen::Vector3d &pt,
        int pt_idx,
        bool pt_is_obs,
        const AABBTree<double,3> *mesh_tree,
        const Eigen::MatrixXd *mesh_x,
        const Eigen::MatrixXi *mesh_tris,
        bool mesh_is_obs,
        double floor_z,
        std::vector<VFCollisionPair> *pairs);
};

// Collision detection against multiple meshes
class EmbeddedMeshCollision : public Collision {
public:
    EmbeddedMeshCollision(const EmbeddedMeshData *mesh_) :
        mesh(mesh_),
        floor_z(-std::numeric_limits<double>::max())
        {}

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
		std::vector<double> *l);

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

/*
class TetMeshCollision : public Collision {
public:
    TetMeshCollision(const TetMeshData *mesh_) :
        mesh(mesh_),
        floor_z(-std::numeric_limits<double>::max())
        {}

    // Performs collision detection.
    // Returns the number of active constraints.
    int detect(
        const Eigen::MatrixXd *x0,
        const Eigen::MatrixXd *x1);

    // Special case for floor since it's common.
    void set_floor(double z)
    {
        floor_z = z;
    }

    // Linearize the constraints and return Jacobian tensor.
    // Constraints are linearized about x for constraint
    // K x = l
    void jacobian(
        const Eigen::MatrixXd *x,
    	std::vector<Eigen::Triplet<double> > *trips_x,
        std::vector<Eigen::Triplet<double> > *trips_y,
    	std::vector<Eigen::Triplet<double> > *trips_z,
		std::vector<double> *l) = 0;

protected:
    const TetMeshData *mesh;
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
*/

} // namespace admmpd

#endif // ADMMPD_COLLISION_H_
