#ifndef UTOPIA_LAPLACIANVIEW_HPP
#define UTOPIA_LAPLACIANVIEW_HPP

#include "utopia_AssemblyView.hpp"

namespace utopia {

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class Laplacian {};


    template<class T>
    class AssemblerView {
    public:
        UTOPIA_INLINE_FUNCTION AssemblerView(const T &value) : value_(value) {}

        template<typename SizeType, class Elem, class Accumulator>
        UTOPIA_INLINE_FUNCTION void assemble(const SizeType &, const Elem &, Accumulator &acc) const {
            acc += value_;
        }

    private:
        const T &value_;
    };

    template<class Mesh, int NComponents, class Quadrature, typename...Args>
    class Laplacian< FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Scalar   = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Elem = typename FunctionSpace::ViewDevice::Elem;
        static const int NNodes = Elem::NNodes;

        using ViewDevice = AssemblerView<StaticMatrix<Scalar, NNodes, NNodes>>;

        Laplacian(const FunctionSpace &space, const Quadrature &q)
        {
            init(space, q);
        }

        ViewDevice view_device() const
        {
            return ViewDevice(mat_);
        }

        template<class Grad, class DX, class Matrix>
        UTOPIA_INLINE_FUNCTION static void assemble(const Grad &grad, const DX &dx, Matrix &mat)
        {
            const auto n = grad.n_points();
            for(SizeType k = 0; k < n; ++k) {
                for(SizeType j = 0; j < grad.n_functions(); ++j) {
                    const auto g_test = grad(j, k);
                    mat(j, j) += dot(g_test, g_test) * dx(k);

                    for(SizeType l = j + 1; l < grad.n_functions(); ++l) {
                        const auto v = dot(g_test, grad(l, k)) * dx(k);
                        mat(j, l) += v;
                        mat(l, j) += v;
                    }
                }
            }
        }

    private:
        StaticMatrix<Scalar, NNodes, NNodes> mat_;

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


#endif //UTOPIA_LAPLACIANVIEW_HPP
