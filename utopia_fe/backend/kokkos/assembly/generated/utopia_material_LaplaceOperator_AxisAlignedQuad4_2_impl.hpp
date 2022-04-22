#ifndef UTOPIA_TPL_MATERIAL_LaplaceOperator_AxisAlignedQuad4_2_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_LaplaceOperator_AxisAlignedQuad4_2_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_AxisAlignedQuad4_2.hpp"
#include "utopia_material_LaplaceOperator.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of LaplaceOperator for symmetric element pair trial=test=AxisAlignedQuad4
         */
        template <typename T, typename GeoT>
        class LaplaceOperator<AxisAlignedQuad4<T, GeoT>> {
        public:
            using ElemT = AxisAlignedQuad4<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "LaplaceOperator<AxisAlignedQuad4>"; }

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
                //	- Result: 4*ADD + 16*ADDAUGMENTEDASSIGNMENT + 12*MUL + 8*POW
                //	- Subexpressions: 8*ADD + 2*DIV + 26*MUL + 4*NEG + 2*POW + 4*SUB
                T x0 = x - 1;
                T x1 = -x0;
                T x2 = py[0] - py[2];
                T x3 = 0.0625 / pow(x2, 2);
                T x4 = y - 1;
                T x5 = -x4;
                T x6 = px[0] - px[2];
                T x7 = 0.0625 / pow(x6, 2);
                T x8 = weight * x2 * x6;
                T x9 = x + 1;
                T x10 = x1 * x3;
                T x11 = x5 * x7;
                T x12 = x8 * (x10 * x9 + x11 * x4);
                T x13 = -x9;
                T x14 = y + 1;
                T x15 = -x14;
                T x16 = x8 * (x10 * x13 + x11 * x15);
                T x17 = x8 * (x0 * x10 + x11 * x14);
                T x18 = x4 * x7;
                T x19 = x3 * x9;
                T x20 = x8 * (x13 * x19 + x15 * x18);
                T x21 = x8 * (x0 * x19 + x14 * x18);
                T x22 = x8 * (x0 * x13 * x3 + x14 * x15 * x7);
                H[0] += x8 * (pow(x1, 2) * x3 + pow(x5, 2) * x7);
                H[1] += x12;
                H[2] += x16;
                H[3] += x17;
                H[4] += x12;
                H[5] += x8 * (x3 * pow(x9, 2) + pow(x4, 2) * x7);
                H[6] += x20;
                H[7] += x21;
                H[8] += x16;
                H[9] += x20;
                H[10] += x8 * (pow(x13, 2) * x3 + pow(x15, 2) * x7);
                H[11] += x22;
                H[12] += x17;
                H[13] += x21;
                H[14] += x22;
                H[15] += x8 * (pow(x0, 2) * x3 + pow(x14, 2) * x7);
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
                // FLOATING POINT OPS!
                //	- Result: 4*ADD + 4*ADDAUGMENTEDASSIGNMENT + 12*MUL
                //	- Subexpressions: 8*ADD + 2*DIV + 19*MUL + 4*NEG + 4*SUB
                T x0 = x + 1;
                T x1 = py[0] - py[2];
                T x2 = 0.25 / x1;
                T x3 = x0 * x2;
                T x4 = x - 1;
                T x5 = x2 * x4;
                T x6 = -x2 * x4;
                T x7 = -x0 * x2;
                T x8 = u[0] * x6 + u[1] * x3 + u[2] * x7 + u[3] * x5;
                T x9 = y - 1;
                T x10 = px[0] - px[2];
                T x11 = 0.25 / x10;
                T x12 = x11 * x9;
                T x13 = y + 1;
                T x14 = x11 * x13;
                T x15 = -x11 * x9;
                T x16 = -x11 * x13;
                T x17 = u[0] * x15 + u[1] * x12 + u[2] * x16 + u[3] * x14;
                T x18 = 4 * weight * x1 * x10;
                Hx[0] += x18 * (x15 * x17 + x6 * x8);
                Hx[1] += x18 * (x12 * x17 + x3 * x8);
                Hx[2] += x18 * (x16 * x17 + x7 * x8);
                Hx[3] += x18 * (x14 * x17 + x5 * x8);
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
                // FLOATING POINT OPS!
                //	- Result: 4*ADD + 4*ADDAUGMENTEDASSIGNMENT + 12*MUL
                //	- Subexpressions: 8*ADD + 2*DIV + 19*MUL + 4*NEG + 4*SUB
                T x0 = x + 1;
                T x1 = py[0] - py[2];
                T x2 = 0.25 / x1;
                T x3 = x0 * x2;
                T x4 = x - 1;
                T x5 = x2 * x4;
                T x6 = -x2 * x4;
                T x7 = -x0 * x2;
                T x8 = u[0] * x6 + u[1] * x3 + u[2] * x7 + u[3] * x5;
                T x9 = y - 1;
                T x10 = px[0] - px[2];
                T x11 = 0.25 / x10;
                T x12 = x11 * x9;
                T x13 = y + 1;
                T x14 = x11 * x13;
                T x15 = -x11 * x9;
                T x16 = -x11 * x13;
                T x17 = u[0] * x15 + u[1] * x12 + u[2] * x16 + u[3] * x14;
                T x18 = 4 * weight * x1 * x10;
                g[0] += x18 * (x15 * x17 + x6 * x8);
                g[1] += x18 * (x12 * x17 + x3 * x8);
                g[2] += x18 * (x16 * x17 + x7 * x8);
                g[3] += x18 * (x14 * x17 + x5 * x8);
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
                //	- Result: 3*ADD + ADDAUGMENTEDASSIGNMENT + 11*MUL + 2*POW
                //	- Subexpressions: 2*ADD + 2*DIV + 4*SUB
                T x0 = px[0] - px[2];
                T x1 = py[0] - py[2];
                T x2 = x + 1;
                T x3 = 1.0 / x1;
                T x4 = x - 1;
                T x5 = y - 1;
                T x6 = 1.0 / x0;
                T x7 = y + 1;
                e += weight * x0 * x1 *
                     (0.0625 * pow(-u[0] * x3 * x4 + u[1] * x2 * x3 - u[2] * x2 * x3 + u[3] * x3 * x4, 2) +
                      0.0625 * pow(-u[0] * x5 * x6 + u[1] * x5 * x6 - u[2] * x6 * x7 + u[3] * x6 * x7, 2));
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
                // FLOATING POINT OPS!
                //	- Result: 11*ADD + 21*ADDAUGMENTEDASSIGNMENT + 27*MUL + 10*POW
                //	- Subexpressions: 14*ADD + 4*DIV + 51*MUL + 4*NEG + 2*POW + 4*SUB
                T x0 = x - 1;
                T x1 = -x0;
                T x2 = py[0] - py[2];
                T x3 = 0.0625 / pow(x2, 2);
                T x4 = y - 1;
                T x5 = -x4;
                T x6 = px[0] - px[2];
                T x7 = 0.0625 / pow(x6, 2);
                T x8 = weight * x2 * x6;
                T x9 = x + 1;
                T x10 = x1 * x3;
                T x11 = x5 * x7;
                T x12 = x8 * (x10 * x9 + x11 * x4);
                T x13 = -x9;
                T x14 = y + 1;
                T x15 = -x14;
                T x16 = x8 * (x10 * x13 + x11 * x15);
                T x17 = x8 * (x0 * x10 + x11 * x14);
                T x18 = x4 * x7;
                T x19 = x3 * x9;
                T x20 = x8 * (x13 * x19 + x15 * x18);
                T x21 = x8 * (x0 * x19 + x14 * x18);
                T x22 = x8 * (x0 * x13 * x3 + x14 * x15 * x7);
                T x23 = 1.0 / x2;
                T x24 = x1 * x23;
                T x25 = x23 * x9;
                T x26 = u[1] * x25;
                T x27 = x0 * x23;
                T x28 = u[3] * x27;
                T x29 = u[0] * x24;
                T x30 = x13 * x23;
                T x31 = u[2] * x30;
                T x32 = 0.0625 * x26 + 0.0625 * x28 + 0.0625 * x29 + 0.0625 * x31;
                T x33 = 1.0 / x6;
                T x34 = x33 * x5;
                T x35 = x33 * x4;
                T x36 = u[1] * x35;
                T x37 = x14 * x33;
                T x38 = u[3] * x37;
                T x39 = u[0] * x34;
                T x40 = x15 * x33;
                T x41 = u[2] * x40;
                T x42 = 0.0625 * x36 + 0.0625 * x38 + 0.0625 * x39 + 0.0625 * x41;
                T x43 = 4 * x8;
                H[0] += x8 * (pow(x1, 2) * x3 + pow(x5, 2) * x7);
                H[1] += x12;
                H[2] += x16;
                H[3] += x17;
                H[4] += x12;
                H[5] += x8 * (x3 * pow(x9, 2) + pow(x4, 2) * x7);
                H[6] += x20;
                H[7] += x21;
                H[8] += x16;
                H[9] += x20;
                H[10] += x8 * (pow(x13, 2) * x3 + pow(x15, 2) * x7);
                H[11] += x22;
                H[12] += x17;
                H[13] += x21;
                H[14] += x22;
                H[15] += x8 * (pow(x0, 2) * x3 + pow(x14, 2) * x7);
                g[0] += x43 * (x24 * x32 + x34 * x42);
                g[1] += x43 * (x25 * x32 + x35 * x42);
                g[2] += x43 * (x30 * x32 + x40 * x42);
                g[3] += x43 * (x27 * x32 + x37 * x42);
                e += x8 * (0.0625 * pow(x26 + x28 + x29 + x31, 2) + 0.0625 * pow(x36 + x38 + x39 + x41, 2));
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using LaplaceOperatorAxisAlignedQuad4 =
            utopia::kokkos::AutoKernel<FE,
                                       utopia::kernels::LaplaceOperator<
                                           utopia::kernels::AxisAlignedQuad4<typename FE::Scalar, typename FE::Scalar>>,
                                       2>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_LaplaceOperator_AxisAlignedQuad4_2_IMPL_hpp
