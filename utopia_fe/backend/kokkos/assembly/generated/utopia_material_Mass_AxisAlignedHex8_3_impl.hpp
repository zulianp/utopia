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
                //	- Subexpressions: 37*MUL + 6*POW + 6*SUB
                T x0 = 1.0 - x;
                T x1 = pow(x0, 2);
                T x2 = 1.0 - y;
                T x3 = pow(x2, 2);
                T x4 = 1.0 - z;
                T x5 = weight * (-px[0] + px[6]) * (-py[0] + py[6]) * (-pz[0] + pz[6]);
                T x6 = pow(x4, 2) * x5;
                T x7 = x3 * x6;
                T x8 = x * x0;
                T x9 = x7 * x8;
                T x10 = x2 * y;
                T x11 = x10 * x6;
                T x12 = x11 * x8;
                T x13 = x1 * x11;
                T x14 = x4 * x5 * z;
                T x15 = x14 * x3;
                T x16 = x1 * x15;
                T x17 = x15 * x8;
                T x18 = x10 * x14;
                T x19 = x18 * x8;
                T x20 = x1 * x18;
                T x21 = pow(x, 2);
                T x22 = x11 * x21;
                T x23 = x15 * x21;
                T x24 = x18 * x21;
                T x25 = pow(y, 2);
                T x26 = x25 * x6;
                T x27 = x26 * x8;
                T x28 = x14 * x25;
                T x29 = x21 * x28;
                T x30 = x28 * x8;
                T x31 = x1 * x28;
                T x32 = x5 * pow(z, 2);
                T x33 = x3 * x32;
                T x34 = x33 * x8;
                T x35 = x10 * x32;
                T x36 = x35 * x8;
                T x37 = x1 * x35;
                T x38 = x21 * x35;
                T x39 = x25 * x32;
                T x40 = x39 * x8;
                H[0] += x1 * x7;
                H[1] += x9;
                H[2] += x12;
                H[3] += x13;
                H[4] += x16;
                H[5] += x17;
                H[6] += x19;
                H[7] += x20;
                H[8] += x9;
                H[9] += x21 * x7;
                H[10] += x22;
                H[11] += x12;
                H[12] += x17;
                H[13] += x23;
                H[14] += x24;
                H[15] += x19;
                H[16] += x12;
                H[17] += x22;
                H[18] += x21 * x26;
                H[19] += x27;
                H[20] += x19;
                H[21] += x24;
                H[22] += x29;
                H[23] += x30;
                H[24] += x13;
                H[25] += x12;
                H[26] += x27;
                H[27] += x1 * x26;
                H[28] += x20;
                H[29] += x19;
                H[30] += x30;
                H[31] += x31;
                H[32] += x16;
                H[33] += x17;
                H[34] += x19;
                H[35] += x20;
                H[36] += x1 * x33;
                H[37] += x34;
                H[38] += x36;
                H[39] += x37;
                H[40] += x17;
                H[41] += x23;
                H[42] += x24;
                H[43] += x19;
                H[44] += x34;
                H[45] += x21 * x33;
                H[46] += x38;
                H[47] += x36;
                H[48] += x19;
                H[49] += x24;
                H[50] += x29;
                H[51] += x30;
                H[52] += x36;
                H[53] += x38;
                H[54] += x21 * x39;
                H[55] += x40;
                H[56] += x20;
                H[57] += x19;
                H[58] += x30;
                H[59] += x31;
                H[60] += x37;
                H[61] += x36;
                H[62] += x40;
                H[63] += x1 * x39;
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
                //	- Subexpressions: 7*ADD + 24*MUL + 6*SUB
                T x0 = 1.0 - x;
                T x1 = 1.0 - y;
                T x2 = 1.0 - z;
                T x3 = x1 * x2;
                T x4 = x0 * x3;
                T x5 = y * z;
                T x6 = x * x5;
                T x7 = x2 * y;
                T x8 = x * x7;
                T x9 = x1 * z;
                T x10 = x * x9;
                T x11 = x0 * x5;
                T x12 = x * x3;
                T x13 = x0 * x7;
                T x14 = x0 * x9;
                T x15 = weight * (-px[0] + px[6]) * (-py[0] + py[6]) * (-pz[0] + pz[6]) *
                        (u[0] * x4 + u[1] * x12 + u[2] * x8 + u[3] * x13 + u[4] * x14 + u[5] * x10 + u[6] * x6 +
                         u[7] * x11);
                Hx[0] += x15 * x4;
                Hx[1] += x12 * x15;
                Hx[2] += x15 * x8;
                Hx[3] += x13 * x15;
                Hx[4] += x14 * x15;
                Hx[5] += x10 * x15;
                Hx[6] += x15 * x6;
                Hx[7] += x11 * x15;
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
                //	- Subexpressions: 7*ADD + 24*MUL + 6*SUB
                T x0 = 1.0 - x;
                T x1 = 1.0 - y;
                T x2 = 1.0 - z;
                T x3 = x1 * x2;
                T x4 = x0 * x3;
                T x5 = y * z;
                T x6 = x * x5;
                T x7 = x2 * y;
                T x8 = x * x7;
                T x9 = x1 * z;
                T x10 = x * x9;
                T x11 = x0 * x5;
                T x12 = x * x3;
                T x13 = x0 * x7;
                T x14 = x0 * x9;
                T x15 = weight * (-px[0] + px[6]) * (-py[0] + py[6]) * (-pz[0] + pz[6]) *
                        (u[0] * x4 + u[1] * x12 + u[2] * x8 + u[3] * x13 + u[4] * x14 + u[5] * x10 + u[6] * x6 +
                         u[7] * x11);
                g[0] += x15 * x4;
                g[1] += x12 * x15;
                g[2] += x15 * x8;
                g[3] += x13 * x15;
                g[4] += x14 * x15;
                g[5] += x10 * x15;
                g[6] += x15 * x6;
                g[7] += x11 * x15;
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
                // TODO
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
                //	- Result: 72*ADDAUGMENTEDASSIGNMENT + 16*MUL
                //	- Subexpressions: 7*ADD + 61*MUL + 6*POW + 6*SUB
                T x0 = 1.0 - x;
                T x1 = pow(x0, 2);
                T x2 = 1.0 - y;
                T x3 = pow(x2, 2);
                T x4 = 1.0 - z;
                T x5 = weight * (-px[0] + px[6]) * (-py[0] + py[6]) * (-pz[0] + pz[6]);
                T x6 = pow(x4, 2) * x5;
                T x7 = x3 * x6;
                T x8 = x * x0;
                T x9 = x7 * x8;
                T x10 = x2 * y;
                T x11 = x10 * x6;
                T x12 = x11 * x8;
                T x13 = x1 * x11;
                T x14 = x1 * x5;
                T x15 = x4 * z;
                T x16 = x15 * x3;
                T x17 = x14 * x16;
                T x18 = x5 * x8;
                T x19 = x16 * x18;
                T x20 = x2 * x4;
                T x21 = x0 * x20;
                T x22 = y * z;
                T x23 = x * x22;
                T x24 = x21 * x23 * x5;
                T x25 = x20 * x22;
                T x26 = x14 * x25;
                T x27 = pow(x, 2);
                T x28 = x11 * x27;
                T x29 = x27 * x5;
                T x30 = x16 * x29;
                T x31 = x25 * x29;
                T x32 = pow(y, 2);
                T x33 = x32 * x6;
                T x34 = x33 * x8;
                T x35 = x15 * x32;
                T x36 = x29 * x35;
                T x37 = x18 * x35;
                T x38 = x14 * x35;
                T x39 = x5 * pow(z, 2);
                T x40 = x3 * x39;
                T x41 = x40 * x8;
                T x42 = x10 * x39;
                T x43 = x42 * x8;
                T x44 = x1 * x42;
                T x45 = x27 * x42;
                T x46 = x32 * x39;
                T x47 = x46 * x8;
                T x48 = x4 * y;
                T x49 = x * x48;
                T x50 = x2 * z;
                T x51 = x * x50;
                T x52 = x0 * x22;
                T x53 = x * x20;
                T x54 = x0 * x48;
                T x55 = x0 * x50;
                T x56 = x5 * (u[0] * x21 + u[1] * x53 + u[2] * x49 + u[3] * x54 + u[4] * x55 + u[5] * x51 + u[6] * x23 +
                              u[7] * x52);
                H[0] += x1 * x7;
                H[1] += x9;
                H[2] += x12;
                H[3] += x13;
                H[4] += x17;
                H[5] += x19;
                H[6] += x24;
                H[7] += x26;
                H[8] += x9;
                H[9] += x27 * x7;
                H[10] += x28;
                H[11] += x12;
                H[12] += x19;
                H[13] += x30;
                H[14] += x31;
                H[15] += x24;
                H[16] += x12;
                H[17] += x28;
                H[18] += x27 * x33;
                H[19] += x34;
                H[20] += x24;
                H[21] += x31;
                H[22] += x36;
                H[23] += x37;
                H[24] += x13;
                H[25] += x12;
                H[26] += x34;
                H[27] += x1 * x33;
                H[28] += x26;
                H[29] += x24;
                H[30] += x37;
                H[31] += x38;
                H[32] += x17;
                H[33] += x19;
                H[34] += x24;
                H[35] += x26;
                H[36] += x1 * x40;
                H[37] += x41;
                H[38] += x43;
                H[39] += x44;
                H[40] += x19;
                H[41] += x30;
                H[42] += x31;
                H[43] += x24;
                H[44] += x41;
                H[45] += x27 * x40;
                H[46] += x45;
                H[47] += x43;
                H[48] += x24;
                H[49] += x31;
                H[50] += x36;
                H[51] += x37;
                H[52] += x43;
                H[53] += x45;
                H[54] += x27 * x46;
                H[55] += x47;
                H[56] += x26;
                H[57] += x24;
                H[58] += x37;
                H[59] += x38;
                H[60] += x44;
                H[61] += x43;
                H[62] += x47;
                H[63] += x1 * x46;
                g[0] += x21 * x56;
                g[1] += x53 * x56;
                g[2] += x49 * x56;
                g[3] += x54 * x56;
                g[4] += x55 * x56;
                g[5] += x51 * x56;
                g[6] += x23 * x56;
                g[7] += x52 * x56;
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
