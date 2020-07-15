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

    void compute_mapping(
        const MatrixXT *V,
        const Eigen::MatrixXi *F);

public:
    SDF();

    bool valid() const { return sdf_dx > 0; };

    // Computes the signed distance field and
    // updates grid mappings
    bool generate(
        const MatrixXT *V,
        const Eigen::MatrixXi *F,
        T dx_frac=-1); // (0-1), fraction of AABB to set dx, -1=auto

    Eigen::AlignedBox<T,3> box() const { return aabb; }
    T dx() const { return sdf_dx; }

    // Samples the SDF at position x in space
    // by mapping it to a cell.
    // < 0 = cell inside surface
    // 0 = cell contains surface
    // > 0 = cell outside surface
    T sample(const VecType &x) const;
    // Sample an area, returns min value
    T sample(const VecType &bmin, const VecType &bmax) const;

    // Given a point on the interior, move along the
    // signed distance field until the surface is reached.
    // Then, project on the surface face.
    // The input mesh (V,F) should be the same as generate(V,F).
    // Returns true if success.
    bool project_out(
        const VecType &pt,
        const MatrixXT *V,
        const Eigen::MatrixXi *F,
        int &face_idx,
        VecType &proj_on_face) const;

    // Returns index into SDF from 3D pt (clamped to AABB bounds)
    Eigen::Vector3i index(const VecType &x) const;

    // Returns list of verts in the cell at ind or box
    void vertices(const Eigen::Vector3i &ind, std::vector<int> &v) const;

    // Returns list of faces in the cell at ind or box
    void faces(const Eigen::Vector3i &ind, std::vector<int> &f) const;
    void faces(const VecType &bmin, const VecType &bmax, std::vector<int> &f) const;

}; // class sdf

} // namespace admmpd

#endif // ADMMPD_SDF_H_
