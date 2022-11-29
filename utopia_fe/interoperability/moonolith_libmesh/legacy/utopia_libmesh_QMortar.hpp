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

#include <libmesh/dense_vector.h>

#include "utopia_libmesh_Utils.hpp"

#include <math.h>
#include <algorithm>
#include <memory>

namespace utopia {

    inline void make_tp(const int i, libMesh::Real &val) {}

    template <class Vec>
    void make_tp(const int i, Vec &val) {
        libMesh::Real s = val(i);
        val.zero();
        val(i) = s;
    }

    void compute_side_normal(const int dim, const libMesh::Elem &side, libMesh::Point &n);

    int order_for_l2_integral(const int dim,
                              const libMesh::Elem &master_el,
                              const int master_order,
                              const libMesh::Elem &slave_el,
                              const int slave_order);

    // void print(const libMesh::QBase &ir, std::ostream &os = std::cout);
    double sum_of_weights(const libMesh::QBase &ir);
    // double sum(const libMesh::DenseMatrix<libMesh::Real> &mat);

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

    bool biorthgonal_weights(const int type, libMesh::Real &w_ii, libMesh::Real &w_ij);

    void mortar_assemble_biorth(const libMesh::FEBase &trial_fe,
                                const libMesh::FEBase &test_fe,
                                const int type,
                                libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble_biorth(const libMesh::FEVectorBase &trial_fe,
                                const libMesh::FEVectorBase &test_fe,
                                const int type,
                                libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble_biorth(const int dim,
                                const libMesh::FEBase &trial_fe,
                                const libMesh::FEBase &test_fe,
                                const int type,
                                const libMesh::DenseVector<libMesh::Real> &indicator,
                                libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble_biorth(const int dim,
                                const libMesh::FEVectorBase &trial_fe,
                                const libMesh::FEVectorBase &test_fe,
                                const int type,
                                const libMesh::DenseVector<libMesh::Real> &indicator,
                                libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble_weights(const libMesh::FEVectorBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights);
    void mortar_assemble_weights(const libMesh::FEBase &fe, libMesh::DenseMatrix<libMesh::Real> &weights);

    void mortar_assemble_weighted_biorth(const libMesh::FEBase &trial_fe,
                                         const libMesh::FEBase &test_fe,
                                         const libMesh::DenseMatrix<libMesh::Real> &weights,
                                         libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble_weighted_biorth(const libMesh::FEVectorBase &trial_fe,
                                         const libMesh::FEVectorBase &test_fe,
                                         const libMesh::DenseMatrix<libMesh::Real> &weights,
                                         libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_normal_and_gap_assemble_weighted_biorth(const libMesh::FEVectorBase &test_fe,
                                                        const int dim,
                                                        const libMesh::Point &surf_normal,
                                                        const libMesh::Point &plane_normal,
                                                        const libMesh::Real &plane_offset,
                                                        const libMesh::DenseMatrix<libMesh::Real> &weights,
                                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                                        libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble_weighted_biorth(const libMesh::FEBase &test_fe,
                                                        const int dim,
                                                        const libMesh::Point &surf_normal,
                                                        const libMesh::Point &plane_normal,
                                                        const libMesh::Real &plane_offset,
                                                        const libMesh::DenseMatrix<libMesh::Real> &weights,
                                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                                        libMesh::DenseVector<libMesh::Real> &gap,
                                                        const bool visdebug = false);

    // void integrate_scalar_function(
    //     const libMesh::FEBase &test_fe,
    //     const std::vector<double> &fun,
    //     libMesh::DenseVector<libMesh::Real> &result
    // );

    template <class Array>
    void integrate_scalar_function(const libMesh::FEBase &test_fe,
                                   const Array &fun,
                                   libMesh::DenseVector<libMesh::Real> &result) {
        const auto &phi = test_fe.get_phi();
        const auto &dx = test_fe.get_JxW();
        const auto n_qp = fun.size();
        const auto n_shape_functions = phi.size();

        assert(n_qp == phi[0].size());
        assert(n_qp == dx.size());

        result.resize(n_shape_functions);
        result.zero();

        for (std::size_t i = 0; i < n_shape_functions; ++i) {
            for (std::size_t qp = 0; qp < n_qp; ++qp) {
                result(i) += phi[i][qp] * fun[qp] * dx[qp];
            }
        }
    }

    void integrate_point_function(const int dim,
                                  const libMesh::FEBase &test_fe,
                                  const std::vector<libMesh::Point> &fun,
                                  libMesh::DenseMatrix<libMesh::Real> &result);

    template <int Dim>
    void l2_project_normal(const libMesh::FEBase &test_fe,
                           const ::moonolith::Vector<double, Dim> &normal,
                           libMesh::DenseVector<libMesh::Real> &result) {
        const auto &phi = test_fe.get_phi();
        const auto &dx = test_fe.get_JxW();
        const auto n_qp = dx.size();
        const auto n_shape_functions = phi.size();

        assert(n_qp == phi[0].size());

        result.resize(n_shape_functions * Dim);
        result.zero();

        for (std::size_t i = 0; i < n_shape_functions; ++i) {
            for (std::size_t qp = 0; qp < n_qp; ++qp) {
                for (std::size_t d = 0; d < Dim; ++d) {
                    result(i * Dim + d) += phi[i][qp] * normal[d] * dx[qp];
                }
            }
        }
    }

    inline double biorth(const std::vector<std::vector<double>> &phi,
                         const libMesh::DenseMatrix<libMesh::Real> &weights,
                         const int i,
                         const int qp) {
        double ret = 0.;
        auto n_funs = phi.size();
        for (std::size_t k = 0; k < n_funs; ++k) {
            ret += phi[k][qp] * weights(i, k);
        }

        return ret;
    }

    template <class Array>
    void integrate_scalar_function_weighted_biorth(const libMesh::FEBase &test_fe,
                                                   const libMesh::DenseMatrix<libMesh::Real> &weights,
                                                   const Array &fun,
                                                   libMesh::DenseVector<libMesh::Real> &result) {
        const auto &phi = test_fe.get_phi();
        const auto &dx = test_fe.get_JxW();
        const auto n_qp = fun.size();
        const auto n_shape_functions = phi.size();

        assert(n_qp == phi[0].size());
        assert(n_qp == dx.size());

        result.resize(n_shape_functions);
        result.zero();

        for (std::size_t i = 0; i < n_shape_functions; ++i) {
            for (std::size_t qp = 0; qp < n_qp; ++qp) {
                auto f = biorth(phi, weights, i, qp);
                result(i) += f * fun[qp] * dx[qp];
            }
        }
    }

    template <int Dim>
    void l2_project_normal_weighted_biorth(const libMesh::FEBase &test_fe,
                                           const libMesh::DenseMatrix<libMesh::Real> &weights,
                                           const ::moonolith::Vector<double, Dim> &normal,

                                           libMesh::DenseVector<libMesh::Real> &result) {
        const auto &phi = test_fe.get_phi();
        const auto &dx = test_fe.get_JxW();
        const auto n_qp = dx.size();
        const auto n_shape_functions = phi.size();

        assert(n_qp == phi[0].size());

        result.resize(n_shape_functions * Dim);
        result.zero();

        for (std::size_t i = 0; i < n_shape_functions; ++i) {
            for (std::size_t qp = 0; qp < n_qp; ++qp) {
                auto f = biorth(phi, weights, i, qp);
                for (std::size_t d = 0; d < Dim; ++d) {
                    result(i * Dim + d) += f * normal[d] * dx[qp];
                }
            }
        }
    }

    void mortar_assemble_weighted_biorth(const libMesh::FEBase &trial_fe,
                                         const libMesh::DenseMatrix<libMesh::Real> &trafo,
                                         const libMesh::FEBase &test_fe,
                                         const libMesh::DenseMatrix<libMesh::Real> &weights,
                                         libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble(const libMesh::FEBase &trial_fe,
                         const libMesh::FEBase &test_fe,
                         libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_assemble(const libMesh::FEVectorBase &trial_fe,
                         const libMesh::FEVectorBase &test_fe,
                         libMesh::DenseMatrix<libMesh::Real> &elmat);

    void mortar_normal_and_gap_assemble(const libMesh::FEBase &test_fe,
                                        const libMesh::DenseVector<libMesh::Real> &surf_normal,
                                        const libMesh::DenseVector<libMesh::Real> &plane_normal,
                                        const libMesh::Real &plane_offset,
                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                        libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble(const uint dim,
                                        const libMesh::FEBase &test_fe,
                                        const libMesh::Point &surf_normal,
                                        const libMesh::Point &plane_normal,
                                        const libMesh::Real &plane_offset,
                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                        libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble(const libMesh::FEVectorBase &test_fe,
                                        const libMesh::DenseVector<libMesh::Real> &surf_normal,
                                        const libMesh::DenseVector<libMesh::Real> &plane_normal,
                                        const libMesh::Real &plane_offset,
                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                        libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble(const uint dim,
                                        const libMesh::FEVectorBase &test_fe,
                                        const libMesh::Point &surf_normal,
                                        const libMesh::Point &plane_normal,
                                        const libMesh::Real &plane_offset,
                                        libMesh::DenseMatrix<libMesh::Real> &normals,
                                        libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble_biorth(const int type,
                                               const libMesh::FEBase &test_fe,
                                               const libMesh::DenseVector<libMesh::Real> &surf_normal,
                                               const libMesh::DenseVector<libMesh::Real> &plane_normal,
                                               const libMesh::Real &plane_offset,
                                               const libMesh::DenseVector<libMesh::Real> &indicator,
                                               libMesh::DenseMatrix<libMesh::Real> &normals,
                                               libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble_biorth(const int type,
                                               const uint dim,
                                               const libMesh::FEBase &test_fe,
                                               const libMesh::Point &surf_normal,
                                               const libMesh::Point &plane_normal,
                                               const libMesh::Real &plane_offset,
                                               const libMesh::DenseVector<libMesh::Real> &indicator,
                                               libMesh::DenseMatrix<libMesh::Real> &normals,
                                               libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble_biorth(const int type,
                                               const libMesh::FEVectorBase &test_fe,
                                               const libMesh::DenseVector<libMesh::Real> &surf_normal,
                                               const libMesh::DenseVector<libMesh::Real> &plane_normal,
                                               const libMesh::Real &plane_offset,
                                               const libMesh::DenseVector<libMesh::Real> &indicator,
                                               libMesh::DenseMatrix<libMesh::Real> &normals,
                                               libMesh::DenseVector<libMesh::Real> &gap);

    void mortar_normal_and_gap_assemble_biorth(const int type,
                                               const uint dim,
                                               const libMesh::FEVectorBase &test_fe,
                                               const libMesh::Point &surf_normal,
                                               const libMesh::Point &plane_normal,
                                               const libMesh::Real &plane_offset,
                                               const libMesh::DenseVector<libMesh::Real> &indicator,
                                               libMesh::DenseMatrix<libMesh::Real> &normals,
                                               libMesh::DenseVector<libMesh::Real> &gap);

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_QMORTAR_HPP
