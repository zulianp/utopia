#ifndef UTOPIA_KOKKOS_WEAK_THERMO_ELASTICITY_COUPLING_HPP
#define UTOPIA_KOKKOS_WEAK_THERMO_ELASTICITY_COUPLING_HPP

#include "utopia_kokkos_LinearElasticityOp.hpp"
#include "utopia_kokkos_Op.hpp"

namespace utopia {
    namespace kokkos {
        namespace kernels {

            template <int Dim, typename Scalar, class Grad, class Function, class Measure>
            class WeakLinearThermoElasticityCouplingOp : public TestTrialOp {
            public:
                using StrainKernel = utopia::kokkos::kernels::LinearizedStrain<Dim, Scalar>;

                UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell,
                                                         const int i,
                                                         const int j,
                                                         const int sub_i) const {
                    Scalar ret = 0.0;
                    for (int qp = 0; qp < n_quad_points; ++qp) {
                        ret += (*this)(cell, i, j, qp, sub_i);
                    }
                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Scalar
                operator()(const int cell, const int i, const int j, const int qp, const int sub_i) const {
                    auto v = -alpha * (3 * lambda + 2 * mu) * fun(j, qp);
                    return StrainKernel::trace(grad, cell, i, qp, sub_i) * v * measure(cell, qp);
                }

                UTOPIA_INLINE_FUNCTION WeakLinearThermoElasticityCouplingOp(const Scalar &lambda,
                                                                            const Scalar &mu,
                                                                            const Scalar &alpha,
                                                                            const Grad &grad,
                                                                            const Function &fun,
                                                                            const Measure &measure)
                    : lambda(lambda),
                      mu(mu),
                      alpha(alpha),
                      grad(grad),
                      fun(fun),
                      measure(measure),
                      n_quad_points(measure.extent(1)) {}

                Scalar lambda, mu, alpha;
                Grad grad;
                Function fun;
                Measure measure;
                int n_quad_points;
            };
        }  // namespace kernels
    }      // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_WEAK_THERMO_ELASTICITY_COUPLING_HPP
