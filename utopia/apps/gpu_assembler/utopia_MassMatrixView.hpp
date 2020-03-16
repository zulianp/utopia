#ifndef UTOPIA_MASSMATRIX_HPP
#define UTOPIA_MASSMATRIX_HPP

#include "utopia_AssemblyView.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_NodalInterpolateView.hpp"
#include "utopia_DeviceOperations.hpp"

namespace utopia {

    template<class Elem>
    class MassMatrixAssembler {
    public:
        using FunValue  = typename Elem::FunValue;
        using Point     = typename Elem::Point;
        using Scalar    = typename Elem::Scalar;

        static const int NFunctions = Elem::NFunctions;

        UTOPIA_INLINE_FUNCTION MassMatrixAssembler(
            const Elem &elem
        )
        : elem_(elem)
        {}

        template<class Quadrature, class Matrix>
        UTOPIA_INLINE_FUNCTION void assemble(const Quadrature &q, Matrix &mat)
        {
            Scalar w;
            Point p;
            for(int qp = 0; qp < Quadrature::NPoints; ++qp) {
                q.point(qp, p);
                w = q.weight(qp) * elem_.measure();

                for(int i = 0; i < NFunctions; ++i) {
                    FunValue &&fun_i = elem_.fun(i, p);
                    mat(i, i) += fun_i * fun_i * w;

                    for(int j = i + 1; j < NFunctions; ++j) {
                        const Scalar v = inner(fun_i, elem_.fun(j, p)) * w;
                        mat(i, j) += v;
                        mat(j, i) += v;
                    }
                }
            }
        }

        template<class Quadrature, class Vector>
        UTOPIA_INLINE_FUNCTION void assemble_vector(const Quadrature &q, Vector &vec)
        {
            Scalar w;
            // Point p;
            for(int qp = 0; qp < Quadrature::NPoints; ++qp) {
                auto && p = q.point(qp);
                w = q.weight(qp) * elem_.measure();

                for(int i = 0; i < NFunctions; ++i) {
                    FunValue &&fun_i = elem_.fun(i, p);
                    vec(i) += fun_i * fun_i * w;

                    for(int j = i + 1; j < NFunctions; ++j) {
                        const Scalar v = inner(fun_i, elem_.fun(j, p)) * w;
                        vec(i) += v;
                        vec(j) += v;
                    }
                }
            }
        }

    private:
        const Elem &elem_;

    };

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class MassMatrix {};

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class ScaledMassMatrix {};

    template<class Mesh, int NComponents, class Quadrature, typename...Args>
    class MassMatrix< FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Scalar   = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Elem = typename FunctionSpace::ViewDevice::Elem;
        static const int NFunctions = Elem::NFunctions;

        using ViewDevice = AssemblerView<StaticMatrix<Scalar, NFunctions, NFunctions>>;

        MassMatrix(const FunctionSpace &space, const Quadrature &q) //: q_(q)
        {
            init(space, q);
        }

        ViewDevice view_device() const
        {
            return ViewDevice(mat_);
        }

        template<class Fun, class DX, class Matrix>
        UTOPIA_INLINE_FUNCTION static void assemble(const Fun &fun, const DX &dx, Matrix &mat)
        {
            const auto n = fun.n_points();
            for(SizeType k = 0; k < n; ++k) {
                for(SizeType j = 0; j < fun.n_functions(); ++j) {
                    const auto g_test = fun(j, k);
                    mat(j, j) += inner(g_test, g_test) * dx(k);

                    for(SizeType l = j + 1; l < fun.n_functions(); ++l) {
                        const auto v = inner(g_test, fun(l, k)) * dx(k);
                        mat(j, l) += v;
                        mat(l, j) += v;
                    }
                }
            }
        }

        void describe(std::ostream &os = std::cout) const
        {
            mat_.describe(os);
        }

    private:
        StaticMatrix<Scalar, NFunctions, NFunctions> mat_;

        void init(const FunctionSpace &space, const Quadrature &q)
        {
            ShapeFunction<FunctionSpace, Quadrature> fun(space, q);
            Differential<FunctionSpace, Quadrature> differential(space, q);

            auto fun_view = fun.view_host();
            auto dx_view  = differential.view_host();

            Elem e;
            space.elem(0, e);

            auto f  = fun_view.make(e);
            auto dx = dx_view.make(e);

            mat_.set(0.0);
            assemble(f, dx, mat_);
        }
    };

    template<class FunctionSpace, class Quadrature>
    MassMatrix<FunctionSpace, Quadrature> mass_matrix(const FunctionSpace &space, const Quadrature &q)
    {
        return MassMatrix<FunctionSpace, Quadrature>(space, q);
    }

    template<class Mesh, int NComponents, class Quadrature, typename...Args>
    class ScaledMassMatrix< FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Vector   = typename FunctionSpace::Vector;
        using Scalar   = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Elem = typename FunctionSpace::ViewDevice::Elem;
        static const int NFunctions = Elem::NFunctions;

        using Differential  = utopia::Differential<FunctionSpace, Quadrature>;
        using Interpolate   = utopia::NodalInterpolate<FunctionSpace, Quadrature>;

        class ViewDevice {
        public:
            using DifferentialView = typename Differential::ViewDevice;
            using InterpolateView  = typename Interpolate::ViewDevice;

            template<typename SizeType, class Elem, class Accumulator>
            UTOPIA_INLINE_FUNCTION void assemble(const SizeType &i, const Elem &e, Accumulator &acc) const {
                auto dx          = dx_.make(i, e);
                auto interpolate = interpolate_.make(i, e);
                auto fun         = interpolate_.fun().make(i, e);

                assemble_aux(fun, dx, interpolate, acc);
            }

            template<typename SizeType, class Elem, class Function, class Accumulator>
            UTOPIA_INLINE_FUNCTION void assemble(const SizeType &i, const Elem &e, Function f, Accumulator &acc) const {
                auto dx          = dx_.make(e);
                auto interpolate = interpolate_.make(e);
                auto fun         = interpolate_.fun().make(e);

                assemble_aux(fun, dx, interpolate, f, acc);
            }

            ViewDevice(const DifferentialView &dx, const InterpolateView &interpolate)
            : dx_(dx), interpolate_(interpolate)
            {}

            DifferentialView dx_;
            InterpolateView interpolate_;
        };

        ScaledMassMatrix(const FunctionSpace &space, const Quadrature &q)
        : dx_(space, q), interpolate_(space, q)
        {}

        void update(const Vector &x)
        {
            interpolate_.update(x);
        }

        ViewDevice view_device() const
        {
            return ViewDevice(dx_.view_device(), interpolate_.view_device());
        }

    private:
        Differential dx_;
        Interpolate interpolate_;

        template<class Fun, class DX, class QPValue, class Function, class Matrix>
        UTOPIA_INLINE_FUNCTION static void assemble_aux(
            const Fun &fun,
            const DX &dx,
            const QPValue &qp_fun,
            Function f,
            Matrix &mat
        )
        {
            const auto n = fun.n_points();
            for(SizeType k = 0; k < n; ++k) {
                for(SizeType j = 0; j < fun.n_functions(); ++j) {
                    const auto g_test = fun(j, k);
                    mat(j, j) += (g_test * g_test) * dx(k) * f(qp_fun(k));

                    for(SizeType l = j + 1; l < fun.n_functions(); ++l) {
                        const auto v = (g_test * fun(l, k)) * dx(k) * f(qp_fun(k));
                        mat(j, l) += v;
                        mat(l, j) += v;
                    }
                }
            }
        }
    };
}



#endif //UTOPIA_MASSMATRIX_HPP
