#ifndef UTOPIA_TPL_MATERIAL_Mass_Line2_1_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_Mass_Line2_1_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Line2_1.hpp"
#include "utopia_material_Mass.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Mass for symmetric element pair trial=test=Line2
         */
        template <typename T, typename GeoT>
        class Mass<Line2<T, GeoT>> {
        public:
            using ElemT = Line2<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Mass<Line2>"; }

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    // TODO
                }

                // TODO
            };

            Mass(const Params &params = Params()) {
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
                //	- Result: 2*ADD + 4*ADDAUGMENTEDASSIGNMENT + 3*MUL + 2*POW
                //	- Subexpressions: ADD + 4*DIV + 3*MUL + 2*SUB
                T x0 = weight * (-1.0 / 2.0 * px[0] + (1.0 / 2.0) * px[1]);
                T x1 = (1.0 / 4.0) * x0;
                T x2 = (1.0 / 2.0) * x;
                T x3 = x0 * (0.5 - x2) * (x2 + 0.5);
                H[0] += x1 * pow(1 - x, 2);
                H[1] += x3;
                H[2] += x3;
                H[3] += x1 * pow(x + 1, 2);
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
                // FLOATING POINT OPS!
                //	- Result: 2*ADDAUGMENTEDASSIGNMENT + 2*MUL
                //	- Subexpressions: 2*ADD + 3*DIV + 4*MUL + 2*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = x0 + 0.5;
                T x3 = weight * (-1.0 / 2.0 * px[0] + (1.0 / 2.0) * px[1]) * (u[0] * x1 + u[1] * x2);
                Hx[0] += x1 * x3;
                Hx[1] += x2 * x3;
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
                // FLOATING POINT OPS!
                //	- Result: 2*ADDAUGMENTEDASSIGNMENT + 2*MUL
                //	- Subexpressions: 2*ADD + 3*DIV + 4*MUL + 2*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = x0 + 0.5;
                T x3 = weight * (-1.0 / 2.0 * px[0] + (1.0 / 2.0) * px[1]) * (u[0] * x1 + u[1] * x2);
                g[0] += x1 * x3;
                g[1] += x2 * x3;
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
                //	- Result: 4*ADD + ADDAUGMENTEDASSIGNMENT + 6*MUL + POW
                //	- Subexpressions: DIV
                T x0 = (1.0 / 2.0) * x;
                e +=
                    weight * (-1.0 / 2.0 * px[0] + (1.0 / 2.0) * px[1]) * pow(u[0] * (0.5 - x0) + u[1] * (x0 + 0.5), 2);
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
                // FLOATING POINT OPS!
                //	- Result: 2*ADD + 7*ADDAUGMENTEDASSIGNMENT + 6*MUL + 3*POW
                //	- Subexpressions: 2*ADD + 4*DIV + 5*MUL + 2*SUB
                T x0 = weight * (-1.0 / 2.0 * px[0] + (1.0 / 2.0) * px[1]);
                T x1 = (1.0 / 4.0) * x0;
                T x2 = (1.0 / 2.0) * x;
                T x3 = x2 + 0.5;
                T x4 = 0.5 - x2;
                T x5 = x0 * x4;
                T x6 = x3 * x5;
                T x7 = u[0] * x4 + u[1] * x3;
                H[0] += x1 * pow(1 - x, 2);
                H[1] += x6;
                H[2] += x6;
                H[3] += x1 * pow(x + 1, 2);
                g[0] += x5 * x7;
                g[1] += x0 * x3 * x7;
                e += x0 * pow(x7, 2);
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using MassLine2 = utopia::kokkos::
            AutoKernel<FE, utopia::kernels::Mass<utopia::kernels::Line2<typename FE::Scalar, typename FE::Scalar>>, 1>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_Mass_Line2_1_IMPL_hpp
