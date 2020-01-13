#ifndef UTOPIA_LINEAR_ELASTICITY_VIEW_HPP
#define UTOPIA_LINEAR_ELASTICITY_VIEW_HPP

#include "utopia_LaplacianView.hpp"

namespace utopia {

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class LinearElasticity {};

    template<class Mesh, int NComponents, class Quadrature, typename...Args>
    class LinearElasticity< FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Scalar   = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Elem = typename FunctionSpace::ViewDevice::Elem;
        static const int NDofs = FunctionSpace::NDofs;

        using ViewDevice = AssemblerView<StaticMatrix<Scalar, NDofs, NDofs>>;

        LinearElasticity(const FunctionSpace &space, const Quadrature &q)
        : mu_(1.0), lambda_(1.0), rescaling_(1.0)
        {
            init(space, q);
        }

        ViewDevice view_device() const
        {
            return ViewDevice(mat_);
        }

        template<class Grad, class DX>
        UTOPIA_INLINE_FUNCTION Scalar kernel(
            const Grad &grad,
            const DX &dx,
            const SizeType &i,
            const SizeType &j,
            const SizeType &qp) const
        {
            auto g_i = grad(i, qp);
            auto g_j = grad(j, qp);

            auto e_i = 0.5 * (transpose(g_i) + g_i);
            auto e_j = 0.5 * (transpose(g_j) + g_j);

            return (
                ((2. * rescaling_) * mu_) * inner(e_i, e_j) +
                (rescaling_ * lambda_) * inner(trace(e_i), trace(e_j))
                ) * dx(qp);
        }

        template<class Grad, class DX, class Matrix>
        UTOPIA_INLINE_FUNCTION void assemble(const Grad &grad, const DX &dx, Matrix &mat) const
        {
            const auto n = grad.n_points();
            for(SizeType k = 0; k < n; ++k) {
                for(SizeType j = 0; j < grad.n_functions(); ++j) {
                    const auto g_test = grad(j, k);
                    mat(j, j) += kernel(grad, dx, j, j, k);

                    for(SizeType l = j + 1; l < grad.n_functions(); ++l) {
                        const auto v = kernel(grad, dx, l, j, k);
                        mat(j, l) += v;
                        mat(l, j) += v;
                    }
                }
            }
        }

    private:
        StaticMatrix<Scalar, NDofs, NDofs> mat_;
        Scalar mu_, lambda_, rescaling_;

        void init(const FunctionSpace &space, const Quadrature &q)
        {
            PhysicalGradient<FunctionSpace, Quadrature> grad(space, q);
            Differential<FunctionSpace, Quadrature> differential(space, q);

            auto grad_view = grad.view_host();
            auto dx_view = differential.view_host();

            Elem e;
            space.elem(0, e);

            auto g  = grad_view.make(0, e);
            auto dx = dx_view.make(0, e);

            mat_.set(0.0);
            assemble(g, dx, mat_);
        }
    };
}

#endif //UTOPIA_LINEAR_ELASTICITY_VIEW_HPP