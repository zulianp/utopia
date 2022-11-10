#ifndef UTOPIA_MORTAR_ASSEMBLE_HPP
#define UTOPIA_MORTAR_ASSEMBLE_HPP

#include <libmesh/elem.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/quadrature_gauss.h>
#include "libmesh/enum_quadrature_type.h"

// // Define the Finite Element object.
#include "libmesh/fe.h"

// // Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

#include "utopia_libmesh_Utils.hpp"

#include "utopia_libmesh_QMortar.hpp"

#include <math.h>
#include <algorithm>
#include <memory>

namespace utopia {

    void compute_side_normal(const int dim, const libMesh::Elem &side, libMesh::Point &n);

    int order_for_l2_integral(const int dim,
                              const libMesh::Elem &master_el,
                              const int master_order,
                              const libMesh::Elem &slave_el,
                              const int slave_order);

    double ref_volume(int type);
    double ref_area_of_surf(int type);

    // class QMortar : public libMesh::QBase {
    // public:
    //     void resize(const int n_points) {
    //         this->get_points().resize(n_points);
    //         this->get_weights().resize(n_points);
    //     }

    //     QMortar(const unsigned int dim, const libMesh::Order order = libMesh::INVALID_ORDER)
    //         : libMesh::QBase(dim, order) {}

    //     libMesh::QuadratureType type() const override { return libMesh::QGAUSS; }

    //     void init_1D(const libMesh::ElemType type = libMesh::INVALID_ELEM, unsigned int p_level = 0) override {
    //         // assert(false);
    //     }

    //     void init_2D(const libMesh::ElemType type = libMesh::INVALID_ELEM, unsigned int p_level = 0) override {
    //         // assert(false);
    //     }

    //     void init_3D(const libMesh::ElemType type = libMesh::INVALID_ELEM, unsigned int p_level = 0) override {
    //         // assert(false);
    //     }
    // };

}  // namespace utopia

#endif  // UTOPIA_MORTAR_ASSEMBLE_HPP
