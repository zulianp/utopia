#ifndef UTOPIA_KOKKOS_NEOHOOKEAN_OP_HPP
#define UTOPIA_KOKKOS_NEOHOOKEAN_OP_HPP

#include "utopia_Base.hpp"
#include "utopia_kokkos_StrainOp.hpp"

namespace utopia {
    namespace kokkos {
        namespace kernels {

            template <int Dim,
                      typename Scalar,
                      class FirstLameParameter,
                      class ShearModulus,
                      class DerformationGradient,
                      class Gradient,
                      class Measure>
            class NeoHookeanOp {
            public:
                using V = StaticVector<Scalar, Dim>;
                using M = StaticMatrix<Scalar, Dim, Dim>;

                // FGF_i = F_inv_t_qp * transpose(Grad_i) * F_inv_t_qp
                UTOPIA_INLINE_FUNCTION static void selective_triple_product(const int di,
                                                                            const M &F_inv_t_qp,
                                                                            const V &grad,
                                                                            M &FGF_i) {
                    // 2x slower
                    // M grad_mat;
                    // make_grad(di, grad, grad_mat);
                    // FGF_i = F_inv_t_qp * transpose(grad_mat) * F_inv_t_qp;

                    M GtF;
                    for (int k1 = 0; k1 < Dim; ++k1) {
                        for (int k2 = 0; k2 < Dim; ++k2) {
                            GtF(k1, k2) = grad[k1] * F_inv_t_qp(di, k2);
                        }
                    }

                    FGF_i = F_inv_t_qp * GtF;
                }

                UTOPIA_INLINE_FUNCTION static Scalar selective_inner(const int di, const M &mat, const V &grad) {
                    Scalar ret = 0.0;

                    for (int dj = 0; dj < Dim; ++dj) {
                        ret += mat(di, dj) * grad[dj];
                    }

                    // 2x slower
                    // M grad_mat;
                    // make_grad(di, grad, grad_mat);
                    // ret = inner(mat, grad_mat);
                    return ret;
                }

                UTOPIA_INLINE_FUNCTION static void make_grad(const int di, const V &vec_grad, M &mat_grad) {
                    mat_grad.set(0.0);

                    for (int i = 0; i < Dim; ++i) {
                        mat_grad(di, i) = vec_grad[i];
                    }
                }

                UTOPIA_INLINE_FUNCTION static void copy_def_grad(const int cell,
                                                                 const int qp,
                                                                 const DerformationGradient &F,
                                                                 M &F_qp) {
                    for (int d1 = 0; d1 < Dim; ++d1) {
                        for (int d2 = 0; d2 < Dim; ++d2) {
                            F_qp(d1, d2) = F(cell, qp, d1 * Dim + d2);
                        }
                    }
                }

                UTOPIA_INLINE_FUNCTION static void linearized_stress(const int sub_trial,
                                                                     const Scalar &lambda,
                                                                     const Scalar &mu,
                                                                     const Scalar &rescale,
                                                                     const V &grad_trial,
                                                                     const M &F_inv_t,
                                                                     const Scalar &J,
                                                                     M &stress_lin) {
                    M FGF;
                    assert(J > 0.0);
                    assert(J != 0.0);
                    assert(J == J);

                    Scalar log_J = device::log(J);
                    assert(log_J == log_J);

                    Scalar beta = (rescale * (lambda * log_J - mu));

                    // FGF = F_inv_t * transpose(Grad_i) * F_inv_t
                    selective_triple_product(sub_trial, F_inv_t, grad_trial, FGF);
                    Scalar inner_F_grad_trial = selective_inner(sub_trial, F_inv_t, grad_trial);

                    make_grad(sub_trial, grad_trial, stress_lin);
                    stress_lin *= rescale * mu;
                    stress_lin += (-beta) * FGF + (lambda * inner_F_grad_trial) * F_inv_t;
                }

                class ApplyHessian {
                public:
                    ApplyHessian(const FirstLameParameter &lambda,
                                 const ShearModulus &mu,
                                 const Scalar &rescale,
                                 const DerformationGradient &F,
                                 const Gradient &grad,
                                 const Measure &measure)
                        : lambda(lambda),
                          mu(mu),
                          rescale(rescale),
                          F(F),
                          grad(grad),
                          measure(measure),
                          n_quad_points(measure.extent(1)) {}

                    UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

                    UTOPIA_INLINE_FUNCTION Scalar
                    operator()(const int cell, const int i, const int j, const int sub_i, const int sub_j) const {
                        Scalar ret = 0.0;
                        for (int qp = 0; qp < n_quad_points; ++qp) {
                            ret += (*this)(cell, i, j, qp, sub_i, sub_j);
                        }

                        return ret;
                    }

                    UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell,
                                                             const int i,
                                                             const int j,
                                                             const int qp,
                                                             const int sub_i,
                                                             const int sub_j) const {
                        V grad_i, grad_j;
                        M stress_lin;
                        M F_qp, F_inv_t_qp;

                        Scalar dX = measure(cell, qp);

                        for (int d = 0; d < Dim; ++d) {
                            grad_i[d] = grad(cell, i, qp, d);
                            grad_j[d] = grad(cell, j, qp, d);
                        }

                        copy_def_grad(cell, qp, F, F_qp);

                        F_inv_t_qp = inv(transpose(F_qp));
                        Scalar J = det(F_qp);
                        linearized_stress(sub_i, lambda, mu, rescale, grad_i, F_inv_t_qp, J, stress_lin);
                        Scalar val = selective_inner(sub_j, stress_lin, grad_j) * dX;
                        return val;
                    }

                    const FirstLameParameter lambda;
                    const ShearModulus mu;
                    const Scalar rescale;
                    const DerformationGradient F;
                    const Gradient grad;
                    const Measure measure;
                    const int n_quad_points;
                };

                template <class Output>
                class StoreHessian {
                public:
                    StoreHessian(const FirstLameParameter &lambda,
                                 const ShearModulus &mu,
                                 const Scalar &rescale,
                                 const DerformationGradient &F,
                                 const Gradient &grad,
                                 const Measure &measure,
                                 Output &data)
                        : lambda(lambda),
                          mu(mu),
                          rescale(rescale),
                          F(F),
                          grad(grad),
                          measure(measure),
                          data(data),
                          n_quad_points(measure.extent(1)) {}

                    UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int i, const int j) const {
                        for (int qp = 0; qp < n_quad_points; ++qp) {
                            (*this)(cell, i, j, qp);
                        }
                    }

                    UTOPIA_INLINE_FUNCTION void operator()(const int cell,
                                                           const int i,
                                                           const int j,
                                                           const int qp) const {
                        V grad_i, grad_j;
                        M FGF_i, stress_lin_i;
                        M F_qp, F_inv_t_qp;

                        auto dX = measure(cell, qp);

                        for (int d = 0; d < Dim; ++d) {
                            grad_i[d] = grad(cell, i, qp, d);
                            grad_j[d] = grad(cell, j, qp, d);
                        }

                        copy_def_grad(cell, qp, F, F_qp);

                        F_inv_t_qp = inv(transpose(F_qp));
                        Scalar J = det(F_qp);
                        assert(J > 0.0);
                        assert(J != 0.0);
                        assert(J == J);

                        Scalar log_J = device::log(J);
                        assert(log_J == log_J);

                        Scalar beta = (rescale * (lambda * log_J - mu));

                        for (int di = 0; di < Dim; ++di) {
                            auto dof_i = i * Dim + di;

                            // FGF_i = F_inv_t_qp * transpose(Grad_i) * F_inv_t_qp
                            selective_triple_product(di, F_inv_t_qp, grad_i, FGF_i);
                            Scalar inner_F_grad_i = selective_inner(di, F_inv_t_qp, grad_i);

                            make_grad(di, grad_i, stress_lin_i);
                            stress_lin_i *= rescale * mu;
                            stress_lin_i += (-beta) * FGF_i + (lambda * inner_F_grad_i) * F_inv_t_qp;

                            for (int dj = 0; dj < Dim; ++dj) {
                                auto dof_j = j * Dim + dj;
                                auto val = selective_inner(dj, stress_lin_i, grad_j) * dX;
                                data(cell, dof_i, dof_j) += val;
                            }
                        }
                    }

                    const FirstLameParameter lambda;
                    const ShearModulus mu;
                    const Scalar rescale;
                    const DerformationGradient F;
                    const Gradient grad;
                    const Measure measure;
                    Output data;
                    const int n_quad_points;
                };

                template <class Output>
                class StoreGradient {
                public:
                    StoreGradient(const FirstLameParameter &lambda,
                                  const ShearModulus &mu,
                                  const Scalar &rescale,
                                  const DerformationGradient &F,
                                  const Gradient &grad,
                                  const Measure &measure,
                                  Output &data)
                        : lambda(lambda),
                          mu(mu),
                          rescale(rescale),
                          F(F),
                          grad(grad),
                          measure(measure),
                          data(data),
                          n_quad_points(measure.extent(1)) {}

                    UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int j) const {
                        for (int qp = 0; qp < n_quad_points; ++qp) {
                            (*this)(cell, j, qp);
                        }
                    }

                    UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int j, const int qp) const {
                        V grad_j;
                        M P;  // P = F_qp
                        M F_inv_t_qp;

                        auto dX = measure(cell, qp);

                        for (int d = 0; d < Dim; ++d) {
                            grad_j[d] = grad(cell, j, qp, d);
                        }

                        copy_def_grad(cell, qp, F, P);

                        F_inv_t_qp = inv(transpose(P));
                        Scalar J = det(P);
                        Scalar log_J = device::log(J);
                        Scalar lambda_log_J = (rescale * lambda) * log_J;

                        assert(J > 0.0);
                        assert(J == J);
                        assert(log_J == log_J);

                        P -= F_inv_t_qp;
                        P *= rescale * mu;
                        P += lambda_log_J * F_inv_t_qp;

                        for (int dj = 0; dj < Dim; ++dj) {
                            auto dof_j = j * Dim + dj;
                            auto val = selective_inner(dj, P, grad_j) * dX;
                            data(cell, dof_j) += val;
                        }
                    }

                    const FirstLameParameter lambda;
                    const ShearModulus mu;
                    const Scalar rescale;
                    const DerformationGradient F;
                    const Gradient grad;
                    const Measure measure;
                    Output data;
                    const int n_quad_points;
                };

                static ApplyHessian apply_hessian_op(const FirstLameParameter &lambda,
                                                     const ShearModulus &mu,
                                                     const Scalar &rescale,
                                                     const DerformationGradient &F,
                                                     const Gradient &grad,
                                                     const Measure &measure) {
                    return ApplyHessian(lambda, mu, rescale, F, grad, measure);
                }

                template <class Output>
                static StoreHessian<Output> store_hessian_op(const FirstLameParameter &lambda,
                                                             const ShearModulus &mu,
                                                             const Scalar &rescale,
                                                             const DerformationGradient &F,
                                                             const Gradient &grad,
                                                             const Measure &measure,
                                                             Output &data) {
                    StoreHessian<Output>(lambda, mu, rescale, F, grad, measure, data);
                }

                template <class Output>
                static StoreGradient<Output> store_gradient_op(const FirstLameParameter &lambda,
                                                               const ShearModulus &mu,
                                                               const Scalar &rescale,
                                                               const DerformationGradient &F,
                                                               const Gradient &grad,
                                                               const Measure &measure,
                                                               Output &data) {
                    StoreGradient<Output>(lambda, mu, rescale, F, grad, measure, data);
                }
            };

        }  // namespace kernels
    }      // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_NEOHOOKEAN_OP_HPP
