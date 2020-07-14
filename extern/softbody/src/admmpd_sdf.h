// Copyright Matt Overby 2020.
// Distributed under the MIT License.

#ifndef ADMMPD_SDF_H_
#define ADMMPD_SDF_H_

#include <Eigen/Geometry>
#include "sdfgen/array3.h"
#include <vector>
#include <unordered_map>
#include <set>

namespace admmpd {

template<typename T>
class SDF {
protected:
    typedef Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> MatrixXT;
    typedef Eigen::Matrix<T,3,1> VecType;

    T sdf_dx;
	sdfgen::Array3f sdf;
    Eigen::AlignedBox<T,3> aabb;
    std::unordered_map<int,std::set<int> > V_map;
    std::unordered_map<int,std::set<int> > F_map;

    // Get a list of cells that overlap with a box
    void get_cells(
        const VecType &bmin, const VecType &bmax,
	    std::vector<int> &cells_inds) const;

public:
    SDF();

    bool valid() const { return sdf_dx > 0; };

    // Computes the signed distance field and
    // updates grid mappings
    bool generate(
        const MatrixXT *V,
        const Eigen::MatrixXi *F,
        T dx_frac=0.1); // (0-1), fraction of AABB to set dx

    // Recomputes mapping of verts and faces to grid
    void update_mapping(
        const MatrixXT *V,
        const Eigen::MatrixXi *F);

    Eigen::AlignedBox<T,3> box() const { return aabb; }

    // Samples the SDF at position x in space
    T sample(const VecType &x) const;

    // Returns index into SDF from 3D pt (clamped to AABB bounds)
    Eigen::Vector3i index(const VecType &x) const;

    // Returns list of verts in the cell at ind
    void vertices(const Eigen::Vector3i &ind, std::vector<int> &v) const;

    // Returns list of faces in the cell at ind
    void faces(const Eigen::Vector3i &ind, std::vector<int> &f) const;

}; // class sdf

} // namespace admmpd

#endif // ADMMPD_SDF_H_
