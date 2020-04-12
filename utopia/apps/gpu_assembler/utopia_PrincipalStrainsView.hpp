#ifndef UTOPIA_PRINCIPAL_STRAINS_VIEW_HPP
#define UTOPIA_PRINCIPAL_STRAINS_VIEW_HPP

#include "utopia_LaplacianView.hpp"
#include "utopia_Utils.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_Coefficient.hpp"

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
        // using SizeType  = typename Elem::SizeType;
        using GradValue = typename Elem::GradValue;
        static const std::size_t NQuadPoints = GradInterpolateView::NQuadPoints;

        PrincipalStrainsView(const GradInterpolateView &grad)
        : grad_(grad)
        {}

        class Evaluation {
        public:
            ArrayView<StaticVector<Scalar, Dim>, NQuadPoints> values;
            ArrayView<GradValue, NQuadPoints> vectors;
            ArrayView<GradValue, NQuadPoints> strain;
        };

        UTOPIA_INLINE_FUNCTION Evaluation make(const Elem &elem) const
        {
            Evaluation ret;

            grad_.get(elem, ret.strain);

            for(std::size_t qp = 0; qp < NQuadPoints; ++qp) {
                ret.strain[qp].symmetrize();
                eig(ret.strain[qp], ret.values[qp], ret.vectors[qp]);
            }

            return ret;
        }

        template<class Matrix>
        UTOPIA_INLINE_FUNCTION static void split_positive(const Evaluation &el_strain, const std::size_t &qp, Matrix &positive)
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
        UTOPIA_INLINE_FUNCTION static void split(const Evaluation &el_strain, const std::size_t &qp, Matrix &negative, Matrix &positive)
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


    private:

        GradInterpolateView grad_;

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


        PrincipalStrains(const std::shared_ptr<Coefficient<FunctionSpace>> &coeff, const Quadrature &q)
        : grad_(coeff, q)
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
