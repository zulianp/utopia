#ifndef UTOPIA_TPL_MATERIAL_LaplaceOperator_Line2_1_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_LaplaceOperator_Line2_1_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Line2_1.hpp"
#include "utopia_material_LaplaceOperator.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of LaplaceOperator for symmetric element pair trial=test=Line2
         */
        template <typename T, typename GeoT>
        class LaplaceOperator<Line2<T, GeoT>> {
        public:
            using ElemT = Line2<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "LaplaceOperator<Line2>"; }

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    // TODO
                }

                // TODO
            };

            LaplaceOperator(const Params &params = Params()) {
                // TODO
            }

            UTOPIA_FUNCTION void hessian(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T weight,
                T *UTOPIA_RESTRICT H) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 4*ADDAUGMENTEDASSIGNMENT
                //	- Subexpressions: DIV + 2*MUL + NEG + SUB
                T x0 = (1.0 / 4.0) * weight / (-1.0 / 2.0 * px[0] + (1.0 / 2.0) * px[1]);
                T x1 = -x0;
                H[0] += x0;
                H[1] += x1;
                H[2] += x1;
                H[3] += x0;
            }

            UTOPIA_FUNCTION void apply(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T weight,
                T *UTOPIA_RESTRICT Hx) const {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                // FLOATING POINT OPS!
                //	- Result: 2*ADDAUGMENTEDASSIGNMENT + MUL
                //	- Subexpressions: 5*DIV + 3*MUL + 2*SUB
                T x0 = 1.0 / (-1.0 / 2.0 * px[0] + (1.0 / 2.0) * px[1]);
                T x1 = weight * (-1.0 / 2.0 * u[0] * x0 + (1.0 / 2.0) * u[1] * x0);
                Hx[0] += -x1;
                Hx[1] += x1;
            }

            UTOPIA_FUNCTION void gradient(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T weight,
                T *UTOPIA_RESTRICT g) const {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                // FLOATING POINT OPS!
                //	- Result: 2*ADDAUGMENTEDASSIGNMENT + MUL
                //	- Subexpressions: 5*DIV + 3*MUL + 2*SUB
                T x0 = 1.0 / (-1.0 / 2.0 * px[0] + (1.0 / 2.0) * px[1]);
                T x1 = weight * (-1.0 / 2.0 * u[0] * x0 + (1.0 / 2.0) * u[1] * x0);
                g[0] += -x1;
                g[1] += x1;
            }

            UTOPIA_FUNCTION void value(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T weight,
                T &e) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: ADD + ADDAUGMENTEDASSIGNMENT + 3*MUL + POW
                //	- Subexpressions: 3*DIV + SUB
                T x0 = -1.0 / 2.0 * px[0] + (1.0 / 2.0) * px[1];
                T x1 = 1.0 / x0;
                e += weight * x0 * pow(-1.0 / 2.0 * u[0] * x1 + (1.0 / 2.0) * u[1] * x1, 2);
            }

            UTOPIA_FUNCTION void eval(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T weight,
                T &e,
                T *UTOPIA_RESTRICT g,
                T *UTOPIA_RESTRICT H) const {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                // FLOATING POINT OPS!
                //	- Result: 7*ADDAUGMENTEDASSIGNMENT + 2*MUL + POW
                //	- Subexpressions: 6*DIV + 4*MUL + NEG + 2*SUB
                T x0 = -1.0 / 2.0 * px[0] + (1.0 / 2.0) * px[1];
                T x1 = 1.0 / x0;
                T x2 = (1.0 / 4.0) * weight * x1;
                T x3 = -x2;
                T x4 = -1.0 / 2.0 * u[0] * x1 + (1.0 / 2.0) * u[1] * x1;
                T x5 = weight * x4;
                H[0] += x2;
                H[1] += x3;
                H[2] += x3;
                H[3] += x2;
                g[0] += -x5;
                g[1] += x5;
                e += weight * x0 * pow(x4, 2);
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using LaplaceOperatorLine2 = utopia::kokkos::AutoKernel<
            FE,
            utopia::kernels::LaplaceOperator<utopia::kernels::Line2<typename FE::Scalar, typename FE::Scalar>>,
            1>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_LaplaceOperator_Line2_1_IMPL_hpp
