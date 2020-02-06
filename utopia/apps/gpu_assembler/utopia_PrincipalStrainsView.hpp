#ifndef UTOPIA_PRINCIPAL_STRAINS_VIEW_HPP
#define UTOPIA_PRINCIPAL_STRAINS_VIEW_HPP

#include "utopia_LaplacianView.hpp"
#include "utopia_Utils.hpp"
#include "utopia_GradInterpolate.hpp"

namespace utopia {

    template<
        class FunctionSpaceView,
        class GradInterpolateView,
        class Elem = typename FunctionSpaceView::Elem,
        class MemType = typename Elem::MemType,
        typename...
    >
    class PrincipalStrainsView {
    public:
        static const int Dim = Elem::Dim;

        using Scalar    = typename Elem::Scalar;
        using SizeType  = typename Elem::SizeType;
        using GradValue = typename Elem::GradValue;
        static const std::size_t NQuadPoints = GradInterpolateView::NQuadPoints;

        PrincipalStrainsView(const GradInterpolateView &grad)
        : grad_(grad)
        {}

        class Evaluation {
        public:
            ArrayView<StaticVector<Scalar, Dim>, NQuadPoints> values;
            ArrayView<GradValue, NQuadPoints> vectors;
        };

        UTOPIA_INLINE_FUNCTION Evaluation make(const Elem &elem) const
        {
            Evaluation ret;
            get(elem, ret.values, ret.vectors);
            return ret;
        }

        template<class Matrix>
        UTOPIA_INLINE_FUNCTION static void split_positive(const Evaluation &el_strain, const SizeType &qp, Matrix &positive)
        {
            positive.set(0.0);

            StaticVector<Scalar, Dim> v;
            for(int d = 0; d < Dim; ++d) {
                auto e_val = el_strain.values[qp][d];
                el_strain.vectors[qp].col(d, v);

                auto outer_v = outer(v, v);

                auto eig_p = split_positive(e_val);

                positive += eig_p * outer_v;
            }
        }

        template<class Matrix>
        UTOPIA_INLINE_FUNCTION static void split(const Evaluation &el_strain, const SizeType &qp, Matrix &negative, Matrix &positive)
        {
            negative.set(0.0);
            positive.set(0.0);

            StaticVector<Scalar, Dim> v;

            for(int d = 0; d < Dim; ++d) {
                auto e_val = el_strain.values[qp][d];
                el_strain.vectors[qp].col(d, v);

                auto outer_v = outer(v, v);

                auto eig_p = split_positive(e_val);
                auto eig_n = split_negative(e_val);

                negative += eig_n * outer_v;
                positive += eig_p * outer_v;
            }
        }

        template<class Values, class Vectors>
        UTOPIA_INLINE_FUNCTION void get(
            const Elem &elem,
            Values &values,
            Vectors &vectors) const
        {
            typename GradInterpolateView::Eval strain;

            grad_.get(elem, strain);

            const SizeType n = strain.size();

            for(SizeType qp = 0; qp < n; ++qp) {
                strain[qp].symmetrize();
                eig(strain[qp], values[qp], vectors[qp]);
                UTOPIA_DEVICE_ASSERT(check(values[qp], vectors[qp]));
            }
        }

    private:

        GradInterpolateView grad_;

        template<class Values, class Vectors>
        UTOPIA_INLINE_FUNCTION static bool check(
            const Values &values,
            const Vectors &vectors)
        {
            auto sum_v = sum(vectors);
            auto sum_e = sum(values);
            UTOPIA_DEVICE_ASSERT( !device::isnan(sum_v) );
            UTOPIA_DEVICE_ASSERT( !device::isnan(sum_e) );

            return !device::isnan(sum_v) && !device::isnan(sum_e);
        }

        UTOPIA_INLINE_FUNCTION static constexpr Scalar split_positive(const Scalar &x) {
            return (device::abs(x) + x)/2;
        }

        UTOPIA_INLINE_FUNCTION static constexpr Scalar split_negative(const Scalar &x) {
            return (device::abs(x) - x)/2;
        }
    };

    template<class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class PrincipalStrains {};

    template<class Mesh, int NComponents, class Quadrature, typename...Args>
    class PrincipalStrains< FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
        using FunctionSpace           = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Vector                  = typename FunctionSpace::Vector;
        using GradInterpolate         = utopia::GradInterpolate<FunctionSpace, Quadrature>;

        using FunctionSpaceViewDevice   = typename FunctionSpace::ViewDevice;
        using GradInterpolateViewDevice = typename GradInterpolate::ViewDevice;

        using ViewDevice              = utopia::PrincipalStrainsView<FunctionSpaceViewDevice, GradInterpolateViewDevice>;


        PrincipalStrains(const FunctionSpace &space, const Quadrature &q)
        : grad_(space, q)
        {}

        inline ViewDevice view_device() const
        {
            return ViewDevice(
                grad_.view_device()
            );
        }

        inline void update(const Vector &x)
        {
            grad_.update(x);
        }

    private:
        GradInterpolate grad_;
    };

}


#endif //UTOPIA_PRINCIPAL_STRAINS_VIEW_HPP
