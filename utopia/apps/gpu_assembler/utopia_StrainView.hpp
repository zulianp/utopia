#ifndef UTOPIA_STRAIN_VIEW_HPP
#define UTOPIA_STRAIN_VIEW_HPP
#include "utopia_FunctionSpace.hpp"
#include "utopia_AssemblyView.hpp"

namespace utopia {


    template<
        class Elem,
        class Quadrature,
        class MemType = typename Elem::MemType,
        typename...
    >
    class StrainView {
    public:
        static const int Dim = Elem::Dim;
        using Scalar = typename Elem::Scalar;
      //        using SizeType = typename Elem::SizeType;
        using Point  = typename Elem::Point;
        using GradValue   = typename Elem::GradValue;

        static const std::size_t NFunctions = Elem::NFunctions;
        static const std::size_t NQPoints   = Quadrature::NPoints;

        UTOPIA_INLINE_FUNCTION StrainView(const ArrayView<GradValue, NQPoints, NFunctions> &values)
        : values_(values)
        {}

        UTOPIA_INLINE_FUNCTION const GradValue &operator()(const int fun_num, const int qp_idx) const
        {
           return values_(qp_idx, fun_num);
        }

        UTOPIA_INLINE_FUNCTION std::size_t n_points() const
        {
            return NQPoints;
        }

        UTOPIA_INLINE_FUNCTION constexpr static std::size_t n_functions()
        {
            return NFunctions;
        }

        UTOPIA_INLINE_FUNCTION const StrainView & make(const Elem &elem) const
        {
            return *this;
        }

    private:
        const ArrayView<GradValue, NQPoints, NFunctions> &values_;

    };

    template<class FunctionSpace, class Quadrature, typename...Args>
    class Strain {
    public:
        using Elem = typename FunctionSpace::ViewDevice::Elem;

        using ViewDevice = utopia::StrainView<Elem, typename Quadrature::ViewDevice>;
        using ViewHost   = utopia::StrainView<Elem, typename Quadrature::ViewHost>;

        using GradValue   = typename Elem::GradValue;

        static const std::size_t NFunctions = Elem::NFunctions;
        static const std::size_t NQPoints   = Quadrature::NPoints;

        Strain(const FunctionSpace &space, const Quadrature &q)
        : values_()
        {
            init(space, q);
        }

        ViewDevice view_device() const
        {
            return ViewDevice(values_);
        }

    private:

        ArrayView<GradValue, NQPoints, NFunctions> values_;

        void init(const FunctionSpace &space, const Quadrature &q)
        {
            PhysicalGradient<FunctionSpace, Quadrature> grad(space, q);

            auto grad_view = grad.view_host();

            Elem e;
            space.elem(0, e);

            auto g  = grad_view.make(e);


            const std::size_t n_qp = q.n_points();
            assert(n_qp == NQPoints);

            for(std::size_t qp = 0; qp < n_qp; ++qp) {
                for(std::size_t i = 0; i < g.n_functions(); ++i) {
                    auto g_i = g(i, qp);
                    values_(qp, i) = 0.5 * (g_i + transpose(g_i));
                }
            }
        }
    };

    //////////////////////////////////////////////////////////////////////////////////
}

#endif //UTOPIA_STRAIN_VIEW_HPP
