#ifndef UTOPIA_LAPLACIANVIEW_HPP
#define UTOPIA_LAPLACIANVIEW_HPP

#include "utopia_AssemblyView.hpp"

namespace utopia {

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class Laplacian {};


    template<typename Scalar>
    class LaplacianKernel {
    public:

        template<class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar apply(
            const Scalar &diff_coeff,
            const Grad   &g_trial,
            const Grad   &g_test,
            const Scalar &dx
            )
        {
            return diff_coeff * inner(g_trial, g_test) * dx;
        }

        template<class Grad>
        UTOPIA_INLINE_FUNCTION static Scalar apply(
            const Scalar &diff_coeff,
            const Scalar &c_trial,
            const Grad   &g_trial,
            const Grad   &g_test,
            const Scalar &dx
            )
        {
            return diff_coeff * c_trial * inner(g_trial, g_test) * dx;
        }

    };

    template<class Elem>
    class LaplacianAssembler {
    public:
        using GradValue = typename Elem::GradValue;
        using Point     = typename Elem::Point;
        using Scalar    = typename Elem::Scalar;

        static const int NFunctions = Elem::NFunctions;

        UTOPIA_INLINE_FUNCTION LaplacianAssembler(
            const Elem &elem,
            const Scalar diffusivity = 1.0)
        : elem_(elem), diffusivity_(diffusivity)
        {}

        UTOPIA_INLINE_FUNCTION void init(const Point &p, const Scalar &weight)
        {
            for(int i = 0; i < NFunctions; ++i) {
                elem_.grad(i, p, g_[i]);
            }

            dx_ = weight * elem_.measure();
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const
        {
            return diffusivity_ * inner(g_[i], g_[j]) * dx_;
        }

        template<class Quadrature, class Matrix>
        UTOPIA_INLINE_FUNCTION void assemble(const Quadrature &q, Matrix &mat)
        {
            Scalar w;
            Point p;
            for(int qp = 0; qp < Quadrature::NPoints; ++qp) {
                q.point(qp, p);
                w = q.weight(qp);

                init(p, w);

                for(int i = 0; i < NFunctions; ++i) {
                    mat(i, i) += (*this)(i, i);

                    for(int j = i + 1; j < NFunctions; ++j) {
                        mat(i, j) += (*this)(i, j);
                    }
                }
            }
        }

    private:
        const Elem &elem_;
        Scalar diffusivity_;
        StaticVector<GradValue, NFunctions> g_;
        Scalar dx_;
    };

    template<class T>
    class AssemblerView {
    public:
        UTOPIA_INLINE_FUNCTION AssemblerView(const T &value) : value_(value) {}

        template<typename SizeType, class Elem, class Accumulator>
        UTOPIA_INLINE_FUNCTION void assemble(const SizeType &, const Elem &, Accumulator &acc) const {
            acc += value_;
        }

        void describe()
        {
            disp(value_);
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
        static const int NDofs = FunctionSpace::NDofs;

        using ViewDevice = AssemblerView<StaticMatrix<Scalar, NDofs, NDofs>>;

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
        StaticMatrix<Scalar, NDofs, NDofs> mat_;

        void init(const FunctionSpace &space, const Quadrature &q)
        {
            PhysicalGradient<FunctionSpace, Quadrature> grad(space, q);
            Differential<FunctionSpace, Quadrature> differential(space, q);

            auto grad_view = grad.view_host();
            auto dx_view = differential.view_host();

            Elem e;
            space.elem(0, e);

            auto g  = grad_view.make(e);
            auto dx = dx_view.make(e);

            mat_.set(0.0);
            assemble(g, dx, mat_);
        }
    };

    template<class FunctionSpace, class Quadrature>
    inline Laplacian<FunctionSpace, Quadrature> laplacian(const FunctionSpace &space, const Quadrature &q)
    {
        return Laplacian<FunctionSpace, Quadrature>(space, q);
    }

}


#endif //UTOPIA_LAPLACIANVIEW_HPP
