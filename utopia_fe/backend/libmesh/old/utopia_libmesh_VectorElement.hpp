#ifndef UTOPIA_LIBMESH_VECTOR_ELEMENT_HPP
#define UTOPIA_LIBMESH_VECTOR_ELEMENT_HPP

#include "utopia_libmesh_FunctionSpace.hpp"

namespace utopia {
    class VectorElement {
    public:
        typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
        typedef TraitsT::JacobianType JacobianType;
        typedef TraitsT::FE FE;

        VectorElement() : grad_flag(false) {}

        template <class Matrix>
        void make_dual(const Matrix &weights) {
            assert(false);
            // auto n_fun = grad.size();
            // auto n_qp = grad[0].size();

            // dual_grad.resize(n_fun);

            // for(uint i = 0; i < n_fun; ++i) {
            //     dual_grad[i].resize(n_qp);

            //     for(uint k = 0; k < n_qp; ++k) {
            //         dual_grad[i][k] = grad[0][k] * weights(i, 0);
            //         for(uint j = 1; j < n_fun; ++j) {
            //             dual_grad[i][k] += grad[j][k] * weights(i, j);
            //         }
            //     }
            // }
        }

        void init_grad(std::vector<std::unique_ptr<FE> > &fe_object) {
            const std::size_t n_quad_points = fe_object[start_var]->get_dphi()[0].size();

            uint n_shape_functions = 0;
            for (uint i = 0; i < n_vars; ++i) {
                n_shape_functions += fe_object[start_var + i]->n_shape_functions();
            }

            grad.resize(n_shape_functions);

            for (std::size_t i = 0; i < n_shape_functions; ++i) {
                grad[i].resize(n_quad_points);
                // TensorValue is by default initialized to 0s
            }

            // NEW
            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                std::size_t offset = 0;
                for (std::size_t i = 0; i < n_vars; ++i) {
                    const auto &fe = fe_object[start_var + i];
                    const uint n_shape_i = fe->n_shape_functions();
                    std::size_t elem_dim = fe->get_dim();

                    for (uint j = 0; j < n_shape_i; ++j, offset++) {
                        grad[offset][qp].resize(elem_dim, dim);
                    }
                }
            }

            reinit_grad(fe_object);
        }

        void reinit_grad(std::vector<std::unique_ptr<FE> > &fe_object) {
            auto n_shape_functions = grad.size();
            auto n_quad_points = grad[0].size();

            for (auto &v : grad) {
                for (auto &w : v) {
                    w.set(0.0);
                }
            }

            for (std::size_t qp = 0; qp < n_quad_points; ++qp) {
                std::size_t offset = 0;

                for (std::size_t i = 0; i < n_vars; ++i) {
                    const auto &fe = fe_object[start_var + i];
                    const uint n_shape_i = fe->n_shape_functions();

                    for (uint j = 0; j < n_shape_i; ++j, offset++) {
                        const auto &grad_i = fe->get_dphi()[j][qp];

                        for (uint d = 0; d < dim; ++d) {
                            grad[offset][qp].set(i, d, grad_i(d));
                        }
                    }
                }
            }
        }

        void init(std::vector<std::unique_ptr<FE> > &fe_object) {
            if (grad_flag) {
                init_grad(fe_object);
            }
        }

        void reinit(std::vector<std::unique_ptr<FE> > &fe_object) {
            if (grad_flag) {
                reinit_grad(fe_object);
            }
        }

        bool grad_flag;

        std::size_t dim;
        std::size_t start_var;
        std::size_t n_vars;
        JacobianType grad;
        // JacobianType dual_grad;
    };
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_VECTOR_ELEMENT_HPP
