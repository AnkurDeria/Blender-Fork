// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_COLLISION_H_
#define ADMMPD_COLLISION_H_

#include "admmpd_bvh.h"
#include "admmpd_types.h"

namespace admmpd {

// I'll update this class/structure another day.
// For now let's get something in place to do floor collisions.
// Probably will work better to use uber-collision class for
// all self and obstacle collisions, reducing the amount of
// for-all vertices loops.
class Collision {
public:
//    virtual void detect(
//        int meshnum,
//        const Eigen::MatrixXd *x,
//        const Eigen::MatrixXi *faces) = 0;
//    virtual void jacobian(
//        const Eigen::MatrixXd *x,
//    	std::vector<Eigen::Triplet<double> > *trips_x,
//        std::vector<Eigen::Triplet<double> > *trips_y,
//    	std::vector<Eigen::Triplet<double> > *trips_z,
//		std::vector<double> *l) = 0;
};

// Collision detection against multiple meshes
class EmbeddedMeshCollision : public Collision {
protected:
    // We progressively build a list of vertices and faces with each
    // add_obstacle call, reindexing as needed. Then we build a tree
    // with all of them. Alternatively we could just build separate trees and combine them.
    Eigen::MatrixXd obs_V0, obs_V1;
    Eigen::MatrixXi obs_F;
    std::vector<Eigen::AlignedBox<double,3> > obs_aabbs;
    AABBTree<double,3> obs_tree;

    Eigen::MatrixXd emb_V0, emb_V1; // copy of embedded vertices
    const Eigen::MatrixXd *emb_barys; // barys of the embedded vtx
    const Eigen::VectorXi *vtx_to_tet; // vertex to tet embedding
    const Eigen::MatrixXi *tets; // tets that embed faces

    struct CollisionPair {
        int p; // point
        int q; // face
        Eigen::Vector3d barys; // barycoords of collision
    };

public:
    // I don't really like having to switch up interface style, but we'll
    // do so here to avoid copies that would happen in admmpd_api.
    void set_obstacles(
        const float *v0,
        const float *v1,
        int nv,
        const int *faces,
        int nf);

    // Updates the tetmesh BVH for self collisions
    // TODO
    void update_bvh(
        const EmbeddedMeshData *mesh,
        const Eigen::MatrixXd *x0,
        const Eigen::MatrixXd *x1)
        { (void)(mesh); (void)(x0); (void)(x1); }

    // Given a list of deformable vertices (the lattice)
    // perform collision detection of the surface mesh against
    // obstacles and possibly self.
    void detect(
        const EmbeddedMeshData *mesh,
        const Eigen::MatrixXd *x0,
        const Eigen::MatrixXd *x1){
            
        }

    void jacobian(
        const Eigen::MatrixXd *x,
    	std::vector<Eigen::Triplet<double> > *trips_x,
        std::vector<Eigen::Triplet<double> > *trips_y,
    	std::vector<Eigen::Triplet<double> > *trips_z,
		std::vector<double> *l)
    {
        
    }
};
/*
class FloorCollider : public Collider {
public:
    virtual void detect(
        int meshnum,
        const Eigen::MatrixXd *x,
        const Eigen::MatrixXi *faces);
    void jacobian(
        const Eigen::MatrixXd *x,
    	std::vector<Eigen::Triplet<double> > *trips_x,
        std::vector<Eigen::Triplet<double> > *trips_y,
    	std::vector<Eigen::Triplet<double> > *trips_z,
		std::vector<double> *l);
};
*/
} // namespace admmpd

#endif // ADMMPD_COLLISION_H_
