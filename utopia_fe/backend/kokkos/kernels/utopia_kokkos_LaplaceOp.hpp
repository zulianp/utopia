#ifndef UTOPIA_KOKKOS_LAPLACE_OP_HPP
#define UTOPIA_KOKKOS_LAPLACE_OP_HPP

#include "utopia_Base.hpp"
#include "utopia_kokkos_Op.hpp"
#include "utopia_kokkos_SubdomainValue.hpp"

namespace utopia {
    namespace kokkos {
        namespace kernels {

            template <typename Scalar, class Coeff, class Gradient, class Measure>
            class LaplaceOp : public TestTrialOp {
            public:
                UTOPIA_INLINE_FUNCTION LaplaceOp(const Coeff &coeff, const Gradient &grad, const Measure &measure)
                    : coeff(coeff),
                      grad(grad),
                      measure(measure),
                      n_qp(measure.extent(1)),
                      spatial_dim_(grad.extent(3)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int i, const int j) const {
                    Scalar ret = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        Scalar dot_g = 0.0;
                        for (int d = 0; d < spatial_dim_; ++d) {
                            dot_g += grad(cell, i, qp, d) * grad(cell, j, qp, d);
                        }

                        ret += dot_g * measure(cell, qp);
                    }

                    return ret * coeff;
                }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int i, const int j, const int qp) const {
                    Scalar ret = 0.0;

                    for (int d = 0; d < spatial_dim_; ++d) {
                        ret += grad(cell, i, qp, d) * grad(cell, j, qp, d);
                    }

                    return ret * coeff * measure(cell, qp);
                }

                UTOPIA_INLINE_FUNCTION static constexpr int dim() { return 1; }

                UTOPIA_INLINE_FUNCTION static constexpr bool is_linear() { return true; }
                UTOPIA_INLINE_FUNCTION static constexpr bool is_operator() { return true; }
                UTOPIA_INLINE_FUNCTION static constexpr bool is_matrix() { return true; }

                UTOPIA_INLINE_FUNCTION static constexpr bool is_vector() { return false; }
                UTOPIA_INLINE_FUNCTION static constexpr bool is_scalar() { return false; }

                const Coeff coeff;
                const Gradient grad;
                const Measure measure;
                const int n_qp;
                const int spatial_dim_;
            };

            template <int Dim, typename Scalar, class Coeff, class Gradient, class Measure>
            class VectorLaplaceOp : public TestTrialOp {
            public:
                UTOPIA_INLINE_FUNCTION Scalar
                operator()(const int cell, const int i, const int j, const int sub_i, const int sub_j) const {
                    if (sub_i != sub_j) {
                        return 0.0;
                    }

                    Scalar ret = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        Scalar val = 0.0;

                        for (int d = 0; d < dim(); ++d) {
                            val += grad(cell, i, qp, d) * grad(cell, j, qp, d);
                        }

                        ret += coeff * val * measure(cell, qp);
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

                UTOPIA_INLINE_FUNCTION VectorLaplaceOp(const Coeff &coeff, const Gradient &grad, const Measure &measure)
                    : coeff(coeff), grad(grad), measure(measure), n_qp(measure.extent(1)) {}

                const Coeff coeff;
                const Gradient grad;
                const Measure measure;
                const int n_qp;
            };

        }  // namespace kernels

        // Handler
        template <class FE, typename DiffusionCoefficient = typename FE::Scalar>
        class LaplaceOp : public TestTrialOpParams, public Configurable {
        public:
            using UniformKernel = utopia::kokkos::kernels::
                LaplaceOp<typename FE::Scalar, DiffusionCoefficient, typename FE::Gradient, typename FE::Measure>;

            using VaryingKernel = utopia::kokkos::SubdomainValueTimesKernel<FE, UniformKernel>;

            void read(Input &in) override {
                TestTrialOpParams::read(in);

                in.get("coeff", coeff);
                in.get("permeability_function", [this](Input &node) {
                    node.get("function", [this](Input &function_node) { subdomain_value.read(function_node); });
                });

                in.get("diffusion_function", [this](Input &node) {
                    node.get("function", [this](Input &function_node) { subdomain_value.read(function_node); });
                });

                bool verbose = false;
                in.get("verbose", verbose);

                if (verbose) {
                    utopia::out() << "------------------------------\n";
                    utopia::out() << name() << ":\n";
                    utopia::out() << "coeff:\t" << coeff << '\n';

                    if (!subdomain_value.empty) {
                        utopia::out() << "function: SubdomainValue\n";
                    }

                    utopia::out() << "------------------------------\n";
                }
            }

            inline static constexpr int n_vars() { return 1; }
            inline static std::string name() { return "LaplaceOperator"; }

            inline bool is_uniform() const { return subdomain_value.empty; }

            inline UniformKernel uniform_kernel(const FE &fe) const {
                return this->update(UniformKernel(coeff, fe.grad(), fe.measure()));
            }

            inline VaryingKernel varying_kernel(const FE &fe) const {
                assert(!is_uniform());
                return VaryingKernel(subdomain_value, uniform_kernel(fe));
            }

            LaplaceOp(const DiffusionCoefficient &coeff = DiffusionCoefficient(1.0))
                : coeff(coeff), subdomain_value(1.0) {}

            DiffusionCoefficient coeff;
            kokkos::SubdomainValue<FE> subdomain_value;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_LAPLACE_OP_HPP
