#ifndef UTOPIA_COEF_STRAINS_VIEW_HPP
#define UTOPIA_COEF_STRAINS_VIEW_HPP

#include "utopia_Coefficient.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_Utils.hpp"

namespace utopia {

    template <class FunctionSpaceView,
              class GradInterpolateView,
              class Elem = typename FunctionSpaceView::Elem,
              class MemType = typename Elem::MemType,
              typename...>
    class CoefStrainView {
    public:
        static const int Dim = Elem::Dim;

        using Scalar = typename Elem::Scalar;
        using GradValue = typename Elem::GradValue;
        static const std::size_t NQuadPoints = GradInterpolateView::NQuadPoints;

        CoefStrainView(const GradInterpolateView &grad) : grad_(grad) {}

        class Evaluation {
        public:
            ArrayView<GradValue, NQuadPoints> strain;
        };

        UTOPIA_INLINE_FUNCTION Evaluation make(const Elem &elem) const {
            Evaluation ret;

            grad_.get(elem, ret.strain);

            for (std::size_t qp = 0; qp < NQuadPoints; ++qp) {
                ret.strain[qp].symmetrize();
            }

            return ret;
        }

    private:
        GradInterpolateView grad_;
    };

    template <class Elem, class Quadrature, class MemType = typename Elem::MemType, typename...>
    class CoefStrain {};

    template <class Mesh, int NComponents, class Quadrature, typename... Args>
    class CoefStrain<FunctionSpace<Mesh, NComponents, Args...>, Quadrature> {
    public:
        using FunctionSpace = utopia::FunctionSpace<Mesh, NComponents, Args...>;
        using Vector = typename FunctionSpace::Vector;
        using GradInterpolate = utopia::GradInterpolate<FunctionSpace, Quadrature>;

        using FunctionSpaceViewDevice = typename FunctionSpace::ViewDevice;
        using GradInterpolateViewDevice = typename GradInterpolate::ViewDevice;

        using ViewDevice = utopia::CoefStrainView<FunctionSpaceViewDevice, GradInterpolateViewDevice>;

        CoefStrain(const std::shared_ptr<Coefficient<FunctionSpace>> &coeff, const Quadrature &q) : grad_(coeff, q) {}

        inline ViewDevice view_device() const { return ViewDevice(grad_.view_device()); }

        inline void update(const Vector &x) { grad_.update(x); }

    private:
        GradInterpolate grad_;
    };

}  // namespace utopia

#endif  // UTOPIA_COEF_STRAINS_VIEW_HPP
