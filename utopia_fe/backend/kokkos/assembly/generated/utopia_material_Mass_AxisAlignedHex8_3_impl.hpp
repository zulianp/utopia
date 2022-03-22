#ifndef UTOPIA_TPL_MATERIAL_Mass_AxisAlignedHex8_3_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_Mass_AxisAlignedHex8_3_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_AxisAlignedHex8_3.hpp"
#include "utopia_material_Mass.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Mass for symmetric element pair trial=test=AxisAlignedHex8
         */
        template <typename T, typename GeoT>
        class Mass<AxisAlignedHex8<T, GeoT>> {
        public:
            using ElemT = AxisAlignedHex8<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Mass<AxisAlignedHex8>"; }

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
                //	- Result: 64*ADDAUGMENTEDASSIGNMENT + 8*MUL
                //	- Subexpressions: 6*ADD + 6*DIV + 44*MUL + 6*POW + 9*SUB
                T x0 = pow(1 - z, 2);
                T x1 = weight * (-px[0] + px[6]) * (-py[0] + py[6]) * (-pz[0] + pz[6]);
                T x2 = (1.0 / 64.0) * x1;
                T x3 = x0 * x2;
                T x4 = pow(1 - x, 2);
                T x5 = pow(1 - y, 2);
                T x6 = x4 * x5;
                T x7 = (1.0 / 16.0) * x1;
                T x8 = x0 * x7;
                T x9 = (1.0 / 2.0) * x;
                T x10 = (0.5 - x9) * (x9 + 0.5);
                T x11 = x10 * x5;
                T x12 = x11 * x8;
                T x13 = (1.0 / 4.0) * x1;
                T x14 = (1.0 / 2.0) * y;
                T x15 = (0.5 - x14) * (x14 + 0.5);
                T x16 = x10 * x15;
                T x17 = x13 * x16;
                T x18 = x0 * x17;
                T x19 = x15 * x8;
                T x20 = x19 * x4;
                T x21 = (1.0 / 2.0) * z;
                T x22 = (0.5 - x21) * (x21 + 0.5);
                T x23 = x22 * x7;
                T x24 = x23 * x6;
                T x25 = x13 * x22;
                T x26 = x11 * x25;
                T x27 = x1 * x16 * x22;
                T x28 = x15 * x25;
                T x29 = x28 * x4;
                T x30 = pow(x + 1, 2);
                T x31 = x3 * x30;
                T x32 = x19 * x30;
                T x33 = x23 * x30;
                T x34 = x33 * x5;
                T x35 = x28 * x30;
                T x36 = pow(y + 1, 2);
                T x37 = x10 * x36;
                T x38 = x37 * x8;
                T x39 = x33 * x36;
                T x40 = x25 * x37;
                T x41 = x36 * x4;
                T x42 = x23 * x41;
                T x43 = pow(z + 1, 2);
                T x44 = x2 * x43;
                T x45 = x43 * x7;
                T x46 = x11 * x45;
                T x47 = x17 * x43;
                T x48 = x15 * x45;
                T x49 = x4 * x48;
                T x50 = x30 * x44;
                T x51 = x30 * x48;
                T x52 = x37 * x45;
                H[0] += x3 * x6;
                H[1] += x12;
                H[2] += x18;
                H[3] += x20;
                H[4] += x24;
                H[5] += x26;
                H[6] += x27;
                H[7] += x29;
                H[8] += x12;
                H[9] += x31 * x5;
                H[10] += x32;
                H[11] += x18;
                H[12] += x26;
                H[13] += x34;
                H[14] += x35;
                H[15] += x27;
                H[16] += x18;
                H[17] += x32;
                H[18] += x31 * x36;
                H[19] += x38;
                H[20] += x27;
                H[21] += x35;
                H[22] += x39;
                H[23] += x40;
                H[24] += x20;
                H[25] += x18;
                H[26] += x38;
                H[27] += x3 * x41;
                H[28] += x29;
                H[29] += x27;
                H[30] += x40;
                H[31] += x42;
                H[32] += x24;
                H[33] += x26;
                H[34] += x27;
                H[35] += x29;
                H[36] += x44 * x6;
                H[37] += x46;
                H[38] += x47;
                H[39] += x49;
                H[40] += x26;
                H[41] += x34;
                H[42] += x35;
                H[43] += x27;
                H[44] += x46;
                H[45] += x5 * x50;
                H[46] += x51;
                H[47] += x47;
                H[48] += x27;
                H[49] += x35;
                H[50] += x39;
                H[51] += x40;
                H[52] += x47;
                H[53] += x51;
                H[54] += x36 * x50;
                H[55] += x52;
                H[56] += x29;
                H[57] += x27;
                H[58] += x40;
                H[59] += x42;
                H[60] += x49;
                H[61] += x47;
                H[62] += x52;
                H[63] += x41 * x44;
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
                // FLOATING POINT OPS!
                //	- Result: 8*ADDAUGMENTEDASSIGNMENT + 8*MUL
                //	- Subexpressions: 10*ADD + 3*DIV + 24*MUL + 6*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * z;
                T x5 = 0.5 - x4;
                T x6 = x3 * x5;
                T x7 = x1 * x6;
                T x8 = x0 + 0.5;
                T x9 = x6 * x8;
                T x10 = x2 + 0.5;
                T x11 = x10 * x5;
                T x12 = x11 * x8;
                T x13 = x1 * x11;
                T x14 = x4 + 0.5;
                T x15 = x14 * x3;
                T x16 = x1 * x15;
                T x17 = x15 * x8;
                T x18 = x10 * x14;
                T x19 = x18 * x8;
                T x20 = x1 * x18;
                T x21 = weight * (-px[0] + px[6]) * (-py[0] + py[6]) * (-pz[0] + pz[6]) *
                        (u[0] * x7 + u[1] * x9 + u[2] * x12 + u[3] * x13 + u[4] * x16 + u[5] * x17 + u[6] * x19 +
                         u[7] * x20);
                Hx[0] += x21 * x7;
                Hx[1] += x21 * x9;
                Hx[2] += x12 * x21;
                Hx[3] += x13 * x21;
                Hx[4] += x16 * x21;
                Hx[5] += x17 * x21;
                Hx[6] += x19 * x21;
                Hx[7] += x20 * x21;
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
                // FLOATING POINT OPS!
                //	- Result: 8*ADDAUGMENTEDASSIGNMENT + 8*MUL
                //	- Subexpressions: 10*ADD + 3*DIV + 24*MUL + 6*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * z;
                T x5 = 0.5 - x4;
                T x6 = x3 * x5;
                T x7 = x1 * x6;
                T x8 = x0 + 0.5;
                T x9 = x6 * x8;
                T x10 = x2 + 0.5;
                T x11 = x10 * x5;
                T x12 = x11 * x8;
                T x13 = x1 * x11;
                T x14 = x4 + 0.5;
                T x15 = x14 * x3;
                T x16 = x1 * x15;
                T x17 = x15 * x8;
                T x18 = x10 * x14;
                T x19 = x18 * x8;
                T x20 = x1 * x18;
                T x21 = weight * (-px[0] + px[6]) * (-py[0] + py[6]) * (-pz[0] + pz[6]) *
                        (u[0] * x7 + u[1] * x9 + u[2] * x12 + u[3] * x13 + u[4] * x16 + u[5] * x17 + u[6] * x19 +
                         u[7] * x20);
                g[0] += x21 * x7;
                g[1] += x21 * x9;
                g[2] += x12 * x21;
                g[3] += x13 * x21;
                g[4] += x16 * x21;
                g[5] += x17 * x21;
                g[6] += x19 * x21;
                g[7] += x20 * x21;
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
                //	- Result: 4*ADD + ADDAUGMENTEDASSIGNMENT + 12*MUL + POW
                //	- Subexpressions: 3*ADD + 3*DIV + 4*MUL + 3*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * z;
                T x5 = 0.5 - x4;
                T x6 = x3 * x5;
                T x7 = x0 + 0.5;
                T x8 = x2 + 0.5;
                T x9 = x5 * x8;
                T x10 = x4 + 0.5;
                T x11 = x10 * x3;
                T x12 = x10 * x8;
                e += weight * (-px[0] + px[6]) * (-py[0] + py[6]) * (-pz[0] + pz[6]) *
                     pow(u[0] * x1 * x6 + u[1] * x6 * x7 + u[2] * x7 * x9 + u[3] * x1 * x9 + u[4] * x1 * x11 +
                             u[5] * x11 * x7 + u[6] * x12 * x7 + u[7] * x1 * x12,
                         2);
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
                // FLOATING POINT OPS!
                //	- Result: 73*ADDAUGMENTEDASSIGNMENT + 17*MUL + POW
                //	- Subexpressions: 13*ADD + 6*DIV + 66*MUL + 6*POW + 9*SUB
                T x0 = pow(1 - z, 2);
                T x1 = weight * (-px[0] + px[6]) * (-py[0] + py[6]) * (-pz[0] + pz[6]);
                T x2 = (1.0 / 64.0) * x1;
                T x3 = x0 * x2;
                T x4 = pow(1 - x, 2);
                T x5 = pow(1 - y, 2);
                T x6 = x4 * x5;
                T x7 = (1.0 / 16.0) * x1;
                T x8 = x0 * x7;
                T x9 = (1.0 / 2.0) * x;
                T x10 = 0.5 - x9;
                T x11 = x9 + 0.5;
                T x12 = x10 * x11;
                T x13 = x12 * x5;
                T x14 = x13 * x8;
                T x15 = (1.0 / 2.0) * y;
                T x16 = 0.5 - x15;
                T x17 = x15 + 0.5;
                T x18 = x16 * x17;
                T x19 = (1.0 / 4.0) * x1;
                T x20 = x12 * x18 * x19;
                T x21 = x0 * x20;
                T x22 = x18 * x8;
                T x23 = x22 * x4;
                T x24 = (1.0 / 2.0) * z;
                T x25 = 0.5 - x24;
                T x26 = x24 + 0.5;
                T x27 = x25 * x26;
                T x28 = x27 * x7;
                T x29 = x28 * x6;
                T x30 = x19 * x27;
                T x31 = x13 * x30;
                T x32 = x16 * x25;
                T x33 = x10 * x32;
                T x34 = x17 * x26;
                T x35 = x11 * x34;
                T x36 = x1 * x33 * x35;
                T x37 = x19 * x32 * x34;
                T x38 = x37 * x4;
                T x39 = pow(x + 1, 2);
                T x40 = x3 * x39;
                T x41 = x22 * x39;
                T x42 = x28 * x39;
                T x43 = x42 * x5;
                T x44 = x37 * x39;
                T x45 = pow(y + 1, 2);
                T x46 = x12 * x45;
                T x47 = x46 * x8;
                T x48 = x42 * x45;
                T x49 = x30 * x46;
                T x50 = x4 * x45;
                T x51 = x28 * x50;
                T x52 = pow(z + 1, 2);
                T x53 = x2 * x52;
                T x54 = x52 * x7;
                T x55 = x13 * x54;
                T x56 = x20 * x52;
                T x57 = x18 * x54;
                T x58 = x4 * x57;
                T x59 = x39 * x53;
                T x60 = x39 * x57;
                T x61 = x46 * x54;
                T x62 = x11 * x32;
                T x63 = x17 * x25;
                T x64 = x11 * x63;
                T x65 = x10 * x63;
                T x66 = x16 * x26;
                T x67 = x10 * x66;
                T x68 = x11 * x66;
                T x69 = x10 * x34;
                T x70 = u[0] * x33 + u[1] * x62 + u[2] * x64 + u[3] * x65 + u[4] * x67 + u[5] * x68 + u[6] * x35 +
                        u[7] * x69;
                T x71 = x1 * x70;
                H[0] += x3 * x6;
                H[1] += x14;
                H[2] += x21;
                H[3] += x23;
                H[4] += x29;
                H[5] += x31;
                H[6] += x36;
                H[7] += x38;
                H[8] += x14;
                H[9] += x40 * x5;
                H[10] += x41;
                H[11] += x21;
                H[12] += x31;
                H[13] += x43;
                H[14] += x44;
                H[15] += x36;
                H[16] += x21;
                H[17] += x41;
                H[18] += x40 * x45;
                H[19] += x47;
                H[20] += x36;
                H[21] += x44;
                H[22] += x48;
                H[23] += x49;
                H[24] += x23;
                H[25] += x21;
                H[26] += x47;
                H[27] += x3 * x50;
                H[28] += x38;
                H[29] += x36;
                H[30] += x49;
                H[31] += x51;
                H[32] += x29;
                H[33] += x31;
                H[34] += x36;
                H[35] += x38;
                H[36] += x53 * x6;
                H[37] += x55;
                H[38] += x56;
                H[39] += x58;
                H[40] += x31;
                H[41] += x43;
                H[42] += x44;
                H[43] += x36;
                H[44] += x55;
                H[45] += x5 * x59;
                H[46] += x60;
                H[47] += x56;
                H[48] += x36;
                H[49] += x44;
                H[50] += x48;
                H[51] += x49;
                H[52] += x56;
                H[53] += x60;
                H[54] += x45 * x59;
                H[55] += x61;
                H[56] += x38;
                H[57] += x36;
                H[58] += x49;
                H[59] += x51;
                H[60] += x58;
                H[61] += x56;
                H[62] += x61;
                H[63] += x50 * x53;
                g[0] += x33 * x71;
                g[1] += x62 * x71;
                g[2] += x64 * x71;
                g[3] += x65 * x71;
                g[4] += x67 * x71;
                g[5] += x68 * x71;
                g[6] += x35 * x71;
                g[7] += x69 * x71;
                e += x1 * pow(x70, 2);
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using MassAxisAlignedHex8 = utopia::kokkos::AutoKernel<
            FE,
            utopia::kernels::Mass<utopia::kernels::AxisAlignedHex8<typename FE::Scalar, typename FE::Scalar>>,
            3>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_Mass_AxisAlignedHex8_3_IMPL_hpp
