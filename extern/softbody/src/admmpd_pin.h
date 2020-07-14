// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_PIN_H_
#define ADMMPD_PIN_H_

#include "admmpd_bvh.h"
#include "admmpd_types.h"
#include "admmpd_embeddedmesh.h"
#include <unordered_map>
#include <vector>

namespace admmpd {

class Pin {
public:

    virtual ~Pin() {}

    // Clears all pins
    virtual void clear() = 0;

    // Set the pin location (q) and per-axis stiffness (k)
    // Stiffness should be 0 to 1. It can go larger, but
    // the resulting matrix will be poorly conditioned.
    virtual void set_pin(
        int idx,
        const Eigen::Vector3d &q,
        const Eigen::Vector3d &k) = 0;

    // Returns true if the vert is pinned
//    virtual bool has_pin(int idx) const = 0;

    // Creates linearization for constraint:
    // Px=q with stiffnesses baked in
    virtual void linearize(
        const Eigen::MatrixXd *x, // not used yet
    	std::vector<Eigen::Triplet<double> > *trips,
		std::vector<double> *q) = 0;

    // Returns per-axis pin stiffness
//    virtual Eigen::Vector3d get_pin_k(int idx) const = 0;

    // Returns pin location, or zero vector if not set
 //   virtual Eigen::Vector3d get_pin_pos(int idx) const = 0;
};

class EmbeddedMeshPin : public Pin {
public:
    EmbeddedMeshPin(const EmbeddedMesh *mesh_) : mesh(mesh_){}

    // Clears all pins
    void clear();

    // Set the pin location of the embedded vertex
    void set_pin(
        int idx,
        const Eigen::Vector3d &p,
        const Eigen::Vector3d &k);

    // Returns true if the deforming vertex is pinned
//    bool has_pin(int idx) const;

    void linearize(
        const Eigen::MatrixXd *x, // not used yet
    	std::vector<Eigen::Triplet<double> > *trips,
		std::vector<double> *q);

    // Returns per-axis pin stiffness of the deforming vertex (not embedded)
    // or zero if not pinned
    // Baryweights are included.
//    Eigen::Vector3d get_pin_k(int idx) const;

    // Returns pin location of the deforming vertex (not embedded)
    // or zero vector if not set
//    Eigen::Vector3d get_pin_pos(int idx) const;

protected:
    // A ptr to the embedded mesh data
    const EmbeddedMesh *mesh;

    // Pins for embedded vertices:
    std::unordered_map<int,Eigen::Vector3d> pin_k;
    std::unordered_map<int,Eigen::Vector3d> pin_pos;
};

} // namespace admmpd

#endif // ADMMPD_PIN_H_
