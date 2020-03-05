#ifndef UTOPIA_MULTI_VARIATE_ELEMENT_HPP
#define UTOPIA_MULTI_VARIATE_ELEMENT_HPP

#include "utopia_Quadrature.hpp"
#include "utopia_ArrayView.hpp"

namespace utopia {

    template<class Elem, int NVariables>
    class MultiVariateElem {
    public:
        using Scalar         = typename Elem::Scalar;
        // using SizeType       = typename Elem::SizeType;
        using MemType        = typename Elem::MemType;
        static const int Dim = Elem::Dim;
        static const int NNodes = Elem::NNodes;
        static const int NFunctions = NNodes * NVariables;

        using UnivGrad      = utopia::StaticVector<Scalar, Dim>;
        using GradValue     = utopia::StaticMatrix<Scalar, NVariables, Dim>;
        // using NodeIndex     = typename Elem::NodeIndex;
        using FunValue      = utopia::StaticVector<Scalar, Dim>;
        using Point         = typename Elem::Point;

        // static_assert(NVariables == Dim, "Number of variables must be equal to dim for vector elements");

        static UTOPIA_INLINE_FUNCTION constexpr SizeType n_functions()
        {
            return NFunctions;
        }

        template<typename Point>
        UTOPIA_INLINE_FUNCTION auto fun(const int i, const Point &p) const -> FunValue
        {
            const int univ_i = i % NNodes;
            const int dim = i / NNodes;
            FunValue f;
            f.set(0.0);
            f[dim] = univar_elem_.fun(univ_i, p);
            return f;
        }

        template<typename Point>
        UTOPIA_INLINE_FUNCTION void node(const std::size_t &i, Point &p) const
        {
            univar_elem_.node(i, p);
        }

        template<typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void centroid(PhysicalPoint &out) const
        {
            univar_elem_.centroid(out);
        }

        template<typename Point, typename OutGrad>
        UTOPIA_INLINE_FUNCTION void grad(const int i, const Point &p, OutGrad &g) const
        {
            const int univ_i = i % NNodes;

            UnivGrad univ_g;
            univar_elem_.grad(univ_i, p, univ_g);

            const int dim = i / NNodes;

            g.set(0.0);

            for(int d = 0; d < Dim; ++d) {
                g(dim, d) = univ_g[d];
            }
        }

        UTOPIA_INLINE_FUNCTION constexpr static bool is_affine()
        {
            return Elem::is_affine();
        }

        UTOPIA_INLINE_FUNCTION constexpr static Scalar reference_measure()
        {
            return Elem::reference_measure();
        }

        UTOPIA_INLINE_FUNCTION Scalar measure() const
        {
            return univar_elem_.measure();
        }

        template<typename RefPoint, typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void point(const RefPoint &in, PhysicalPoint &out) const
        {
            univar_elem_.point(in, out);
        }

        template<class... Args>
        UTOPIA_INLINE_FUNCTION MultiVariateElem(Args &&...args)
        : univar_elem_(std::forward<Args>(args)...)
        {}

        template<class H>
        UTOPIA_INLINE_FUNCTION void set(
            const StaticVector<Scalar, Dim> &translation,
            const H &h)
        {
           univar_elem_.set(translation, h);
        }

        UTOPIA_INLINE_FUNCTION const StaticVector2<Scalar> &translation() const
        {
            return univar_elem_.translation_();
        }

        Elem &univar_elem()
        {
            return univar_elem_;
        }

        const Elem &univar_elem() const
        {
            return univar_elem_;
        }

        inline SizeType idx() const
        {
            return univar_elem_.idx();
        }

        inline bool is_valid() const
        {
            return univar_elem_.is_valid();
        }

    private:
        Elem univar_elem_;
    };

    template<class Elem>
    class MultiVariateElem<Elem, 1> final : public Elem {
    public:
        static const int NFunctions = Elem::NNodes;

        Elem &univar_elem()
        {
            return *this;
        }

        const Elem &univar_elem() const
        {
            return *this;
        }
    };

    template<typename Elem, int NVariables, int Dim, int Order, typename ...Args>
    class Quadrature<MultiVariateElem<Elem, NVariables>, Order, Dim, Args...> : public Quadrature<Elem, Order, Dim, Args...>
    {};

}

#endif //UTOPIA_MULTI_VARIATE_ELEMENT_HPP
