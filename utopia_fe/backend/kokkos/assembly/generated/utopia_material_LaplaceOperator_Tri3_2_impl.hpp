#ifndef UTOPIA_TPL_MATERIAL_LaplaceOperator_Tri3_2_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_LaplaceOperator_Tri3_2_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Tri3_2.hpp"
#include "utopia_material_LaplaceOperator.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of LaplaceOperator for symmetric element pair trial=test=Tri3
         */
        template <typename T, typename GeoT>
        class LaplaceOperator<Tri3<T, GeoT>> {
        public:
            using ElemT = Tri3<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "LaplaceOperator<Tri3>"; }

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
                const GeoT *UTOPIA_RESTRICT py,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T weight,
                T *UTOPIA_RESTRICT H) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 3*ADD + 9*ADDAUGMENTEDASSIGNMENT + 9*MUL + 6*POW
                //	- Subexpressions: 6*ADD + 18*MUL + 6*NEG + POW + 5*SUB
                T x0 = -px[2];
                T x1 = -px[1] - x0;
                T x2 = px[0] - px[1];
                T x3 = -py[2];
                T x4 = py[0] + x3;
                T x5 = px[0] + x0;
                T x6 = py[0] - py[1];
                T x7 = pow(x2 * x4 - x5 * x6, -2);
                T x8 = py[1] + x3;
                T x9 = -x2;
                T x10 = -x4;
                T x11 = weight * (x10 * x9 - x5 * x6);
                T x12 = x1 * x7;
                T x13 = x7 * x8;
                T x14 = x11 * (x10 * x13 + x12 * x5);
                T x15 = x11 * (x12 * x9 + x13 * x6);
                T x16 = x11 * (x10 * x6 * x7 + x5 * x7 * x9);
                H[0] += x11 * (pow(x1, 2) * x7 + x7 * pow(x8, 2));
                H[1] += x14;
                H[2] += x15;
                H[3] += x14;
                H[4] += x11 * (pow(x10, 2) * x7 + pow(x5, 2) * x7);
                H[5] += x16;
                H[6] += x15;
                H[7] += x16;
                H[8] += x11 * (pow(x6, 2) * x7 + x7 * pow(x9, 2));
            }

            UTOPIA_FUNCTION void apply(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T weight,
                T *UTOPIA_RESTRICT Hx) const {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                // FLOATING POINT OPS!
                //	- Result: 3*ADD + 3*ADDAUGMENTEDASSIGNMENT + 9*MUL
                //	- Subexpressions: 7*ADD + DIV + 18*MUL + 5*NEG + 5*SUB
                T x0 = px[0] - px[1];
                T x1 = -py[2];
                T x2 = py[0] + x1;
                T x3 = -px[2];
                T x4 = px[0] + x3;
                T x5 = py[0] - py[1];
                T x6 = 1.0 / (x0 * x2 - x4 * x5);
                T x7 = x6 * (-px[1] - x3);
                T x8 = u[1] * x6;
                T x9 = -x0;
                T x10 = u[2] * x6;
                T x11 = u[0] * x7 + x10 * x9 + x4 * x8;
                T x12 = x6 * (py[1] + x1);
                T x13 = -x2;
                T x14 = u[0] * x12 + x10 * x5 + x13 * x8;
                T x15 = 3 * weight * (x13 * x9 - x4 * x5);
                T x16 = x11 * x6;
                T x17 = x14 * x6;
                Hx[0] += x15 * (x11 * x7 + x12 * x14);
                Hx[1] += x15 * (x13 * x17 + x16 * x4);
                Hx[2] += x15 * (x16 * x9 + x17 * x5);
            }

            UTOPIA_FUNCTION void gradient(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
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
                UTOPIA_UNUSED(y);
                // FLOATING POINT OPS!
                //	- Result: 3*ADD + 3*ADDAUGMENTEDASSIGNMENT + 9*MUL
                //	- Subexpressions: 7*ADD + DIV + 18*MUL + 5*NEG + 5*SUB
                T x0 = px[0] - px[1];
                T x1 = -py[2];
                T x2 = py[0] + x1;
                T x3 = -px[2];
                T x4 = px[0] + x3;
                T x5 = py[0] - py[1];
                T x6 = 1.0 / (x0 * x2 - x4 * x5);
                T x7 = x6 * (-px[1] - x3);
                T x8 = u[1] * x6;
                T x9 = -x0;
                T x10 = u[2] * x6;
                T x11 = u[0] * x7 + x10 * x9 + x4 * x8;
                T x12 = x6 * (py[1] + x1);
                T x13 = -x2;
                T x14 = u[0] * x12 + x10 * x5 + x13 * x8;
                T x15 = 3 * weight * (x13 * x9 - x4 * x5);
                T x16 = x11 * x6;
                T x17 = x14 * x6;
                g[0] += x15 * (x11 * x7 + x12 * x14);
                g[1] += x15 * (x13 * x17 + x16 * x4);
                g[2] += x15 * (x16 * x9 + x17 * x5);
            }

            UTOPIA_FUNCTION void value(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T weight,
                T &e) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 6*ADD + ADDAUGMENTEDASSIGNMENT + 11*MUL + 2*POW
                //	- Subexpressions: 2*ADD + DIV + 5*MUL + 4*NEG + 3*SUB
                T x0 = px[0] - px[1];
                T x1 = -x0;
                T x2 = -py[2];
                T x3 = py[0] + x2;
                T x4 = -x3;
                T x5 = -px[2];
                T x6 = px[0] + x5;
                T x7 = py[0] - py[1];
                T x8 = 1.0 / (x0 * x3 - x6 * x7);
                T x9 = u[0] * x8;
                T x10 = u[1] * x8;
                T x11 = u[2] * x8;
                e += weight * (x1 * x4 - x6 * x7) *
                     (pow(x1 * x11 + x10 * x6 + x9 * (-px[1] - x5), 2) +
                      pow(x10 * x4 + x11 * x7 + x9 * (py[1] + x2), 2));
            }

            UTOPIA_FUNCTION void eval(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T weight,
                T &e,
                T *UTOPIA_RESTRICT g,
                T *UTOPIA_RESTRICT H) const {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                // FLOATING POINT OPS!
                //	- Result: 7*ADD + 13*ADDAUGMENTEDASSIGNMENT + 19*MUL + 8*POW
                //	- Subexpressions: 10*ADD + DIV + 31*MUL + 6*NEG + POW + 5*SUB
                T x0 = -px[2];
                T x1 = -px[1] - x0;
                T x2 = px[0] - px[1];
                T x3 = -py[2];
                T x4 = py[0] + x3;
                T x5 = px[0] + x0;
                T x6 = py[0] - py[1];
                T x7 = x2 * x4 - x5 * x6;
                T x8 = pow(x7, -2);
                T x9 = py[1] + x3;
                T x10 = -x2;
                T x11 = -x4;
                T x12 = weight * (x10 * x11 - x5 * x6);
                T x13 = x1 * x8;
                T x14 = x8 * x9;
                T x15 = x12 * (x11 * x14 + x13 * x5);
                T x16 = x12 * (x10 * x13 + x14 * x6);
                T x17 = x12 * (x10 * x5 * x8 + x11 * x6 * x8);
                T x18 = 1.0 / x7;
                T x19 = x1 * x18;
                T x20 = u[1] * x18;
                T x21 = u[2] * x18;
                T x22 = u[0] * x19 + x10 * x21 + x20 * x5;
                T x23 = x18 * x9;
                T x24 = u[0] * x23 + x11 * x20 + x21 * x6;
                T x25 = 3 * x12;
                T x26 = x18 * x22;
                T x27 = x18 * x24;
                H[0] += x12 * (pow(x1, 2) * x8 + x8 * pow(x9, 2));
                H[1] += x15;
                H[2] += x16;
                H[3] += x15;
                H[4] += x12 * (pow(x11, 2) * x8 + pow(x5, 2) * x8);
                H[5] += x17;
                H[6] += x16;
                H[7] += x17;
                H[8] += x12 * (pow(x10, 2) * x8 + pow(x6, 2) * x8);
                g[0] += x25 * (x19 * x22 + x23 * x24);
                g[1] += x25 * (x11 * x27 + x26 * x5);
                g[2] += x25 * (x10 * x26 + x27 * x6);
                e += x12 * (pow(x22, 2) + pow(x24, 2));
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using LaplaceOperatorTri3 = utopia::kokkos::AutoKernel<
            FE,
            utopia::kernels::LaplaceOperator<utopia::kernels::Tri3<typename FE::Scalar, typename FE::Scalar>>,
            2>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_LaplaceOperator_Tri3_2_IMPL_hpp
