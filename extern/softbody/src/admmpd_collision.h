// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_COLLISION_H_
#define ADMMPD_COLLISION_H_

#include "admmpd_bvh.h"
#include "admmpd_types.h"
#include "admmpd_embeddedmesh.h"
#include <set>

namespace admmpd {

struct VFCollisionPair {
    int p_idx; // point
    int p_is_obs; // 0 or 1
    int q_idx; // face, or -1 if floor
    int q_is_obs; // 0 or 1
    Eigen::Vector3d q_pt; // pt of collision (if q obs)
    Eigen::Vector3d q_bary; // barys of collision (if q !obs)
    VFCollisionPair();
};

class Collision {
public:

    struct Settings {
        double collision_eps;
        double floor_z;
        bool test_floor;
        bool self_collision;
        Settings() :
            collision_eps(1e-10),
            floor_z(-std::numeric_limits<double>::max()),
            test_floor(true),
            self_collision(false)
            {}
    } settings;

    struct ObstacleData {
        bool has_obs() const { return F.rows()>0; }
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        AABBTree<double,3> tree;
        std::vector<Eigen::AlignedBox<double,3> > leaves;
    } obsdata;

    virtual ~Collision() {}

    // Performs collision detection.
    // Returns the number of active constraints.
    virtual int detect(
        const Eigen::MatrixXd *x0,
        const Eigen::MatrixXd *x1) = 0;

    // Appends the per-vertex graph of dependencies
    // for constraints (ignores obstacles).
    virtual void graph(
        std::vector<std::set<int> > &g) = 0;

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
    virtual void set_floor(double z) { settings.floor_z=z; }
    virtual double get_floor() const { return settings.floor_z; }

    // Linearize the constraints and return Jacobian.
    virtual void linearize(
        const Eigen::MatrixXd *x,
    	std::vector<Eigen::Triplet<double> > *trips,
		std::vector<double> *d) = 0;

    // Given a point and a mesh, perform
    // discrete collision detection.
    std::pair<bool,VFCollisionPair>
    detect_against_obs(
        const Eigen::Vector3d &pt,
        const ObstacleData *obs) const;

    // Given a point and a mesh, perform
    // discrete collision detection.
//    std::pair<bool,VFCollisionPair>
//    detect_point_against_mesh(
//        int pt_idx,
//        bool pt_is_obs,
//        const Eigen::Vector3d &pt,
//        bool mesh_is_obs,
//        const EmbeddedMesh *emb_mesh,
//        const Eigen::MatrixXd *mesh_tets_x,
//        const AABBTree<double,3> *mesh_tets_tree) const;
};

// Collision detection against multiple meshes
class EmbeddedMeshCollision : public Collision {
public:
    EmbeddedMeshCollision(const EmbeddedMesh *mesh_) : mesh(mesh_){}

    // Performs collision detection and stores pairs
    int detect(
        const Eigen::MatrixXd *x0,
        const Eigen::MatrixXd *x1);

    void graph(
        std::vector<std::set<int> > &g);
    
    // Linearizes the collision pairs about x
    // for the constraint Kx=l
    void linearize(
        const Eigen::MatrixXd *x,
    	std::vector<Eigen::Triplet<double> > *trips,
		std::vector<double> *d);

protected:
    // A ptr to the embedded mesh data
    const EmbeddedMesh *mesh;

    // Pairs are compute on detect
    std::vector<Eigen::Vector2i> vf_pairs; // index into per_vertex_pairs
    std::vector<std::vector<VFCollisionPair> > per_vertex_pairs;

    // Updates the tetmesh BVH for self collisions.
    // Called by detect()
    // TODO
    void update_bvh(
        const Eigen::MatrixXd *x0,
        const Eigen::MatrixXd *x1);
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
