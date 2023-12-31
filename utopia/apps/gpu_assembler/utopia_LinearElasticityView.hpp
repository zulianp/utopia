#ifndef UTOPIA_LINEAR_ELASTICITY_VIEW_HPP
#define UTOPIA_LINEAR_ELASTICITY_VIEW_HPP

#include "utopia_AppBase.hpp"
#include "utopia_LaplacianView.hpp"

namespace utopia {

    template <class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class LinearElasticity {};

    template <typename Scalar>
    class LinearElasticityKernel {
    public:
        template <class Strain, class DX>
        UTOPIA_INLINE_FUNCTION static auto strain_apply(const Scalar &mu,
                                                        const Scalar &lambda,
                                                        const Strain &e_i,
                                                        const Strain &e_j,
                                                        const DX &dx) -> decltype(dx * lambda) {
            return ((2. * mu) * inner(e_i, e_j) + lambda * trace(e_i) * trace(e_j)) * dx;
        }

        template <class Grad, class DX>
        UTOPIA_INLINE_FUNCTION static auto apply(const Scalar &mu,
                                                 const Scalar &lambda,
                                                 const Grad &g_i,
                                                 const Grad &g_j,
                                                 const DX &dx) -> decltype(dx * lambda) {
            auto e_i = 0.5 * (transpose(g_i) + g_i);
            auto e_j = 0.5 * (transpose(g_j) + g_j);
            return strain_apply(mu, lambda, e_i, e_j, dx);
        }

        template <class Grad, class DX>
        UTOPIA_INLINE_FUNCTION static auto apply(const Scalar &mu,
                                                 const Scalar &lambda,
                                                 const Scalar &c_i,
                                                 const Grad &g_i,
                                                 const Grad &g_j,
                                                 const DX &dx) -> decltype(dx * lambda) {
            return apply(mu, lambda, g_i, g_j, dx) * c_i;
        }
    };

    template <class Mesh, int NComponents, class Quadrature, typename... Args>
    class LinearElasticity<FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Elem = typename FunctionSpace::ViewDevice::Elem;
        static const int NDofs = FunctionSpace::NDofs;

        using ViewDevice = AssemblerView<StaticMatrix<Scalar, NDofs, NDofs>>;

        LinearElasticity(const FunctionSpace &space,
                         const Quadrature &q,
                         const Scalar &mu = 1.0,
                         const Scalar &lambda = 1.0)
            : mu_(mu), lambda_(lambda), rescaling_(1.0) {
            init(space, q);
        }

        ViewDevice view_device() const { return ViewDevice(mat_); }

        template <class Grad, class DX>
        UTOPIA_INLINE_FUNCTION auto kernel(const Grad &grad,
                                           const DX &dx,
                                           const SizeType &i,
                                           const SizeType &j,
                                           const SizeType &qp) const -> decltype(dx(0)) {
            return LinearElasticityKernel<Scalar>::apply(
                rescaling_ * mu_, rescaling_ * lambda_, grad(i, qp), grad(j, qp), dx(qp));
        }

        template <class Grad, class DX, class Matrix>
        UTOPIA_INLINE_FUNCTION void assemble(const Grad &grad, const DX &dx, Matrix &mat) const {
            const int n = grad.n_points();
            const int n_fun = grad.n_functions();

            for (int k = 0; k < n; ++k) {
                // pragma GCCunroll(NDofs)
                for (int j = 0; j < n_fun; ++j) {
                    // const auto g_test = grad(j, k);
                    mat(j, j) += simd::integrate(kernel(grad, dx, j, j, k));

                    for (int l = j + 1; l < n_fun; ++l) {
                        const auto v = simd::integrate(kernel(grad, dx, l, j, k));
                        mat(j, l) += v;
                        mat(l, j) += v;
                    }
                }
            }
        }

    private:
        StaticMatrix<Scalar, NDofs, NDofs> mat_;
        Scalar mu_, lambda_, rescaling_;

        void init(const FunctionSpace &space, const Quadrature &q) {
            PhysicalGradient<FunctionSpace, Quadrature> grad(space, q);
            Differential<FunctionSpace, Quadrature> differential(space, q);

            auto grad_view = grad.view_host();
            auto dx_view = differential.view_host();

            Elem e;
            space.elem(0, e);

            auto g = grad_view.make(e);
            auto dx = dx_view.make(e);

            mat_.set(0.0);
            assemble(g, dx, mat_);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_LINEAR_ELASTICITY_VIEW_HPP