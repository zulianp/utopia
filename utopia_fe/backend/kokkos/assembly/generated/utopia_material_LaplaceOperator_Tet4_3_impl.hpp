#ifndef UTOPIA_TPL_MATERIAL_LaplaceOperator_Tet4_3_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_LaplaceOperator_Tet4_3_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Tet4_3.hpp"
#include "utopia_material_LaplaceOperator.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of LaplaceOperator for symmetric element pair trial=test=Tet4
         */
        template <typename T, typename GeoT>
        class LaplaceOperator<Tet4<T, GeoT>> {
        public:
            using ElemT = Tet4<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "LaplaceOperator<Tet4>"; }

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
                const GeoT *UTOPIA_RESTRICT pz,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T weight,
                T *UTOPIA_RESTRICT H) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 4*ADD + 16*ADDAUGMENTEDASSIGNMENT + 13*MUL + 12*POW
                //	- Subexpressions: 14*ADD + DIV + 65*MUL + 4*NEG + POW + 27*SUB
                T x0 = -px[0] + px[1];
                T x1 = -py[0] + py[2];
                T x2 = x0 * x1;
                T x3 = -px[0] + px[2];
                T x4 = -py[0] + py[1];
                T x5 = x3 * x4;
                T x6 = x2 - x5;
                T x7 = -pz[0] + pz[3];
                T x8 = -pz[0] + pz[1];
                T x9 = -py[0] + py[3];
                T x10 = x3 * x9;
                T x11 = -pz[0] + pz[2];
                T x12 = -px[0] + px[3];
                T x13 = x0 * x9;
                T x14 = x1 * x12;
                T x15 = x10 * x8 + x11 * x12 * x4 - x11 * x13 - x14 * x8 + x2 * x7 - x5 * x7;
                T x16 = 1.0 / x15;
                T x17 = x16 * x6;
                T x18 = x12 * x4 - x13;
                T x19 = x16 * x18;
                T x20 = x10 - x14;
                T x21 = x16 * x20;
                T x22 = -x17 - x19 - x21;
                T x23 = -x0 * x11 + x3 * x8;
                T x24 = x16 * x23;
                T x25 = x0 * x7 - x12 * x8;
                T x26 = x16 * x25;
                T x27 = x11 * x12 - x3 * x7;
                T x28 = x16 * x27;
                T x29 = -x24 - x26 - x28;
                T x30 = -x1 * x8 + x11 * x4;
                T x31 = x16 * x30;
                T x32 = -x4 * x7 + x8 * x9;
                T x33 = x16 * x32;
                T x34 = x1 * x7 - x11 * x9;
                T x35 = x16 * x34;
                T x36 = -x31 - x33 - x35;
                T x37 = weight * x15;
                T x38 = x37 * (x21 * x22 + x28 * x29 + x35 * x36);
                T x39 = x37 * (x19 * x22 + x26 * x29 + x33 * x36);
                T x40 = x37 * (x17 * x22 + x24 * x29 + x31 * x36);
                T x41 = pow(x15, -2);
                T x42 = x20 * x41;
                T x43 = x27 * x41;
                T x44 = x34 * x41;
                T x45 = x37 * (x18 * x42 + x25 * x43 + x32 * x44);
                T x46 = x37 * (x23 * x43 + x30 * x44 + x42 * x6);
                T x47 = x37 * (x18 * x41 * x6 + x23 * x25 * x41 + x30 * x32 * x41);
                H[0] += x37 * (pow(x22, 2) + pow(x29, 2) + pow(x36, 2));
                H[1] += x38;
                H[2] += x39;
                H[3] += x40;
                H[4] += x38;
                H[5] += x37 * (pow(x20, 2) * x41 + pow(x27, 2) * x41 + pow(x34, 2) * x41);
                H[6] += x45;
                H[7] += x46;
                H[8] += x39;
                H[9] += x45;
                H[10] += x37 * (pow(x18, 2) * x41 + pow(x25, 2) * x41 + pow(x32, 2) * x41);
                H[11] += x47;
                H[12] += x40;
                H[13] += x46;
                H[14] += x47;
                H[15] += x37 * (pow(x23, 2) * x41 + pow(x30, 2) * x41 + x41 * pow(x6, 2));
            }

            UTOPIA_FUNCTION void apply(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T weight,
                T *UTOPIA_RESTRICT Hx) const {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                UTOPIA_UNUSED(z);
                // FLOATING POINT OPS!
                //	- Result: 4*ADD + 4*ADDAUGMENTEDASSIGNMENT + 16*MUL
                //	- Subexpressions: 11*ADD + DIV + 48*MUL + 3*NEG + 27*SUB
                T x0 = -px[0] + px[1];
                T x1 = -py[0] + py[2];
                T x2 = x0 * x1;
                T x3 = -px[0] + px[2];
                T x4 = -py[0] + py[1];
                T x5 = x3 * x4;
                T x6 = -pz[0] + pz[3];
                T x7 = -pz[0] + pz[1];
                T x8 = -py[0] + py[3];
                T x9 = x3 * x8;
                T x10 = -pz[0] + pz[2];
                T x11 = -px[0] + px[3];
                T x12 = x0 * x8;
                T x13 = x1 * x11;
                T x14 = x10 * x11 * x4 - x10 * x12 - x13 * x7 + x2 * x6 - x5 * x6 + x7 * x9;
                T x15 = 1.0 / x14;
                T x16 = x15 * (x2 - x5);
                T x17 = x15 * (x11 * x4 - x12);
                T x18 = x15 * (-x13 + x9);
                T x19 = -x16 - x17 - x18;
                T x20 = u[0] * x19 + u[1] * x18 + u[2] * x17 + u[3] * x16;
                T x21 = x15 * (-x0 * x10 + x3 * x7);
                T x22 = x15 * (x0 * x6 - x11 * x7);
                T x23 = x15 * (x10 * x11 - x3 * x6);
                T x24 = -x21 - x22 - x23;
                T x25 = u[0] * x24 + u[1] * x23 + u[2] * x22 + u[3] * x21;
                T x26 = x15 * (-x1 * x7 + x10 * x4);
                T x27 = x15 * (-x4 * x6 + x7 * x8);
                T x28 = x15 * (x1 * x6 - x10 * x8);
                T x29 = -x26 - x27 - x28;
                T x30 = u[0] * x29 + u[1] * x28 + u[2] * x27 + u[3] * x26;
                T x31 = 4 * weight * x14;
                Hx[0] += x31 * (x19 * x20 + x24 * x25 + x29 * x30);
                Hx[1] += x31 * (x18 * x20 + x23 * x25 + x28 * x30);
                Hx[2] += x31 * (x17 * x20 + x22 * x25 + x27 * x30);
                Hx[3] += x31 * (x16 * x20 + x21 * x25 + x26 * x30);
            }

            UTOPIA_FUNCTION void gradient(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T weight,
                T *UTOPIA_RESTRICT g) const {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                UTOPIA_UNUSED(z);
                // FLOATING POINT OPS!
                //	- Result: 4*ADD + 4*ADDAUGMENTEDASSIGNMENT + 16*MUL
                //	- Subexpressions: 11*ADD + DIV + 48*MUL + 3*NEG + 27*SUB
                T x0 = -px[0] + px[1];
                T x1 = -py[0] + py[2];
                T x2 = x0 * x1;
                T x3 = -px[0] + px[2];
                T x4 = -py[0] + py[1];
                T x5 = x3 * x4;
                T x6 = -pz[0] + pz[3];
                T x7 = -pz[0] + pz[1];
                T x8 = -py[0] + py[3];
                T x9 = x3 * x8;
                T x10 = -pz[0] + pz[2];
                T x11 = -px[0] + px[3];
                T x12 = x0 * x8;
                T x13 = x1 * x11;
                T x14 = x10 * x11 * x4 - x10 * x12 - x13 * x7 + x2 * x6 - x5 * x6 + x7 * x9;
                T x15 = 1.0 / x14;
                T x16 = x15 * (x2 - x5);
                T x17 = x15 * (x11 * x4 - x12);
                T x18 = x15 * (-x13 + x9);
                T x19 = -x16 - x17 - x18;
                T x20 = u[0] * x19 + u[1] * x18 + u[2] * x17 + u[3] * x16;
                T x21 = x15 * (-x0 * x10 + x3 * x7);
                T x22 = x15 * (x0 * x6 - x11 * x7);
                T x23 = x15 * (x10 * x11 - x3 * x6);
                T x24 = -x21 - x22 - x23;
                T x25 = u[0] * x24 + u[1] * x23 + u[2] * x22 + u[3] * x21;
                T x26 = x15 * (-x1 * x7 + x10 * x4);
                T x27 = x15 * (-x4 * x6 + x7 * x8);
                T x28 = x15 * (x1 * x6 - x10 * x8);
                T x29 = -x26 - x27 - x28;
                T x30 = u[0] * x29 + u[1] * x28 + u[2] * x27 + u[3] * x26;
                T x31 = 4 * weight * x14;
                g[0] += x31 * (x19 * x20 + x24 * x25 + x29 * x30);
                g[1] += x31 * (x18 * x20 + x23 * x25 + x28 * x30);
                g[2] += x31 * (x17 * x20 + x22 * x25 + x27 * x30);
                g[3] += x31 * (x16 * x20 + x21 * x25 + x26 * x30);
            }

            UTOPIA_FUNCTION void value(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T weight,
                T &e) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 7*ADD + ADDAUGMENTEDASSIGNMENT + 22*MUL + 3*POW
                //	- Subexpressions: 2*ADD + DIV + 34*MUL + 21*SUB
                T x0 = -pz[0] + pz[3];
                T x1 = -px[0] + px[1];
                T x2 = -py[0] + py[2];
                T x3 = x1 * x2;
                T x4 = -pz[0] + pz[1];
                T x5 = -px[0] + px[2];
                T x6 = -py[0] + py[3];
                T x7 = x5 * x6;
                T x8 = -pz[0] + pz[2];
                T x9 = -px[0] + px[3];
                T x10 = -py[0] + py[1];
                T x11 = x1 * x6;
                T x12 = x10 * x5;
                T x13 = x2 * x9;
                T x14 = -x0 * x12 + x0 * x3 + x10 * x8 * x9 - x11 * x8 - x13 * x4 + x4 * x7;
                T x15 = 1.0 / x14;
                T x16 = x15 * (-x13 + x7);
                T x17 = x15 * (x10 * x9 - x11);
                T x18 = x15 * (-x12 + x3);
                T x19 = x15 * (-x0 * x5 + x8 * x9);
                T x20 = x15 * (x0 * x1 - x4 * x9);
                T x21 = x15 * (-x1 * x8 + x4 * x5);
                T x22 = x15 * (x0 * x2 - x6 * x8);
                T x23 = x15 * (-x0 * x10 + x4 * x6);
                T x24 = x15 * (x10 * x8 - x2 * x4);
                e += weight * x14 *
                     (pow(u[0] * (-x16 - x17 - x18) + u[1] * x16 + u[2] * x17 + u[3] * x18, 2) +
                      pow(u[0] * (-x19 - x20 - x21) + u[1] * x19 + u[2] * x20 + u[3] * x21, 2) +
                      pow(u[0] * (-x22 - x23 - x24) + u[1] * x22 + u[2] * x23 + u[3] * x24, 2));
            }

            UTOPIA_FUNCTION void eval(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T weight,
                T &e,
                T *UTOPIA_RESTRICT g,
                T *UTOPIA_RESTRICT H) const {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                UTOPIA_UNUSED(z);
                // FLOATING POINT OPS!
                //	- Result: 9*ADD + 21*ADDAUGMENTEDASSIGNMENT + 30*MUL + 15*POW
                //	- Subexpressions: 23*ADD + DIV + 78*MUL + 4*NEG + POW + 27*SUB
                T x0 = -px[0] + px[1];
                T x1 = -py[0] + py[2];
                T x2 = x0 * x1;
                T x3 = -px[0] + px[2];
                T x4 = -py[0] + py[1];
                T x5 = x3 * x4;
                T x6 = x2 - x5;
                T x7 = -pz[0] + pz[3];
                T x8 = -pz[0] + pz[1];
                T x9 = -py[0] + py[3];
                T x10 = x3 * x9;
                T x11 = -pz[0] + pz[2];
                T x12 = -px[0] + px[3];
                T x13 = x0 * x9;
                T x14 = x1 * x12;
                T x15 = x10 * x8 + x11 * x12 * x4 - x11 * x13 - x14 * x8 + x2 * x7 - x5 * x7;
                T x16 = 1.0 / x15;
                T x17 = x16 * x6;
                T x18 = x12 * x4 - x13;
                T x19 = x16 * x18;
                T x20 = x10 - x14;
                T x21 = x16 * x20;
                T x22 = -x17 - x19 - x21;
                T x23 = -x0 * x11 + x3 * x8;
                T x24 = x16 * x23;
                T x25 = x0 * x7 - x12 * x8;
                T x26 = x16 * x25;
                T x27 = x11 * x12 - x3 * x7;
                T x28 = x16 * x27;
                T x29 = -x24 - x26 - x28;
                T x30 = -x1 * x8 + x11 * x4;
                T x31 = x16 * x30;
                T x32 = -x4 * x7 + x8 * x9;
                T x33 = x16 * x32;
                T x34 = x1 * x7 - x11 * x9;
                T x35 = x16 * x34;
                T x36 = -x31 - x33 - x35;
                T x37 = weight * x15;
                T x38 = x37 * (x21 * x22 + x28 * x29 + x35 * x36);
                T x39 = x37 * (x19 * x22 + x26 * x29 + x33 * x36);
                T x40 = x37 * (x17 * x22 + x24 * x29 + x31 * x36);
                T x41 = pow(x15, -2);
                T x42 = x20 * x41;
                T x43 = x27 * x41;
                T x44 = x34 * x41;
                T x45 = x37 * (x18 * x42 + x25 * x43 + x32 * x44);
                T x46 = x37 * (x23 * x43 + x30 * x44 + x42 * x6);
                T x47 = x37 * (x18 * x41 * x6 + x23 * x25 * x41 + x30 * x32 * x41);
                T x48 = u[0] * x22 + u[1] * x21 + u[2] * x19 + u[3] * x17;
                T x49 = u[0] * x29 + u[1] * x28 + u[2] * x26 + u[3] * x24;
                T x50 = u[0] * x36 + u[1] * x35 + u[2] * x33 + u[3] * x31;
                T x51 = 4 * x37;
                H[0] += x37 * (pow(x22, 2) + pow(x29, 2) + pow(x36, 2));
                H[1] += x38;
                H[2] += x39;
                H[3] += x40;
                H[4] += x38;
                H[5] += x37 * (pow(x20, 2) * x41 + pow(x27, 2) * x41 + pow(x34, 2) * x41);
                H[6] += x45;
                H[7] += x46;
                H[8] += x39;
                H[9] += x45;
                H[10] += x37 * (pow(x18, 2) * x41 + pow(x25, 2) * x41 + pow(x32, 2) * x41);
                H[11] += x47;
                H[12] += x40;
                H[13] += x46;
                H[14] += x47;
                H[15] += x37 * (pow(x23, 2) * x41 + pow(x30, 2) * x41 + x41 * pow(x6, 2));
                g[0] += x51 * (x22 * x48 + x29 * x49 + x36 * x50);
                g[1] += x51 * (x21 * x48 + x28 * x49 + x35 * x50);
                g[2] += x51 * (x19 * x48 + x26 * x49 + x33 * x50);
                g[3] += x51 * (x17 * x48 + x24 * x49 + x31 * x50);
                e += x37 * (pow(x48, 2) + pow(x49, 2) + pow(x50, 2));
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FunctionSpace, class FE>
        using LaplaceOperatorTet4 = utopia::kokkos::AutoKernel<FunctionSpace,
            FE,
            utopia::kernels::LaplaceOperator<utopia::kernels::Tet4<typename FE::Scalar, typename FE::Scalar>>,
            3>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_LaplaceOperator_Tet4_3_IMPL_hpp
