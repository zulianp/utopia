#ifndef UTOPIA_LIBMESH_QMORTAR_HPP
#define UTOPIA_LIBMESH_QMORTAR_HPP

#include "moonolith_map_quadrature.hpp"

#include <libmesh/elem.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/quadrature_gauss.h>
#include "libmesh/enum_quadrature_type.h"

// // Define the Finite Element object.
#include "libmesh/fe.h"

// // Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

#include "utopia_libmesh_Utils.hpp"

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

    class QMortar : public libMesh::QBase {
    public:
        void resize(const int n_points) {
            this->get_points().resize(n_points);
            this->get_weights().resize(n_points);
        }

        QMortar(const unsigned int dim, const libMesh::Order order = libMesh::INVALID_ORDER)
            : libMesh::QBase(dim, order) {}

        libMesh::QuadratureType type() const override { return libMesh::QGAUSS; }

        void init_1D(const libMesh::ElemType type = libMesh::INVALID_ELEM, unsigned int p_level = 0) override {
            // assert(false);
        }

        void init_2D(const libMesh::ElemType type = libMesh::INVALID_ELEM, unsigned int p_level = 0) override {
            // assert(false);
        }

        void init_3D(const libMesh::ElemType type = libMesh::INVALID_ELEM, unsigned int p_level = 0) override {
            // assert(false);
        }
    };

    template <typename T, int Dim>
    inline void convert(const ::moonolith::Quadrature<T, Dim> &q_in,
                        const ::moonolith::Vector<T, Dim> &point_shift,
                        const T &point_rescale,
                        const T &weight_rescale,
                        QMortar &q_out) {
        const auto &p_in = q_in.points;
        const auto &w_in = q_in.weights;

        auto &p_out = q_out.get_points();
        auto &w_out = q_out.get_weights();

        const std::size_t n_qp = q_in.n_points();

        q_out.resize(n_qp);

        for (std::size_t k = 0; k < n_qp; ++k) {
            w_out[k] = w_in[k] * weight_rescale;

            const auto &pk_in = p_in[k];
            auto &pk_out = p_out[k];

            for (int i = 0; i < Dim; ++i) {
                pk_out(i) = pk_in[i] * point_rescale + point_shift[i];
            }
        }
    }

    template <typename T, int Dim>
    inline void convert(const ::moonolith::Quadrature<T, Dim> &q_in, const T &weight_rescale, QMortar &q_out) {
        static const ::moonolith::Vector<T, Dim> zero;
        return convert(q_in, zero, 1., weight_rescale, q_out);
    }

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_QMORTAR_HPP
