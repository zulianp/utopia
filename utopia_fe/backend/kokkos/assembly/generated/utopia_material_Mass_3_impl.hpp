#ifndef UTOPIA_TPL_MATERIAL_Mass_3_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_Mass_3_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Hex8_3.hpp"
#include "utopia_material_Mass.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Mass for symmetric element pair trial=test=Hex8
         */
        template <typename T, typename GeoT>
        class Mass<Hex8<T, GeoT>> {
        public:
            using ElemT = Hex8<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Mass<Hex8>"; }

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
                //	- Subexpressions: 29*ADD + 173*MUL + 6*POW + 42*SUB
                T x0 = 1.0 - x;
                T x1 = pow(x0, 2);
                T x2 = 1.0 - y;
                T x3 = pow(x2, 2);
                T x4 = 1.0 - z;
                T x5 = x * y;
                T x6 = x * x2;
                T x7 = x0 * y;
                T x8 = x0 * x2;
                T x9 = -px[0] * x8 - px[1] * x6 - px[2] * x5 - px[3] * x7 + px[4] * x0 * x2 + px[5] * x * x2 +
                       px[6] * x * y + px[7] * x0 * y;
                T x10 = x * z;
                T x11 = x * x4;
                T x12 = x0 * z;
                T x13 = x0 * x4;
                T x14 = -py[0] * x13 - py[1] * x11 + py[2] * x * x4 + py[3] * x0 * x4 - py[4] * x12 - py[5] * x10 +
                        py[6] * x * z + py[7] * x0 * z;
                T x15 = y * z;
                T x16 = x4 * y;
                T x17 = x2 * z;
                T x18 = x2 * x4;
                T x19 = -pz[0] * x18 + pz[1] * x2 * x4 + pz[2] * x4 * y - pz[3] * x16 - pz[4] * x17 + pz[5] * x2 * z +
                        pz[6] * y * z - pz[7] * x15;
                T x20 = -px[0] * x13 - px[1] * x11 + px[2] * x * x4 + px[3] * x0 * x4 - px[4] * x12 - px[5] * x10 +
                        px[6] * x * z + px[7] * x0 * z;
                T x21 = -py[0] * x18 + py[1] * x2 * x4 + py[2] * x4 * y - py[3] * x16 - py[4] * x17 + py[5] * x2 * z +
                        py[6] * y * z - py[7] * x15;
                T x22 = -pz[0] * x8 - pz[1] * x6 - pz[2] * x5 - pz[3] * x7 + pz[4] * x0 * x2 + pz[5] * x * x2 +
                        pz[6] * x * y + pz[7] * x0 * y;
                T x23 = -px[0] * x18 + px[1] * x2 * x4 + px[2] * x4 * y - px[3] * x16 - px[4] * x17 + px[5] * x2 * z +
                        px[6] * y * z - px[7] * x15;
                T x24 = -py[0] * x8 - py[1] * x6 - py[2] * x5 - py[3] * x7 + py[4] * x0 * x2 + py[5] * x * x2 +
                        py[6] * x * y + py[7] * x0 * y;
                T x25 = -pz[0] * x13 - pz[1] * x11 + pz[2] * x * x4 + pz[3] * x0 * x4 - pz[4] * x12 - pz[5] * x10 +
                        pz[6] * x * z + pz[7] * x0 * z;
                T x26 = weight * (-x14 * x19 * x9 + x14 * x22 * x23 + x19 * x20 * x24 - x20 * x21 * x22 +
                                  x21 * x25 * x9 - x23 * x24 * x25);
                T x27 = x26 * pow(x4, 2);
                T x28 = x27 * x3;
                T x29 = x * x0;
                T x30 = x28 * x29;
                T x31 = x5 * x8;
                T x32 = x27 * x31;
                T x33 = x2 * y;
                T x34 = x27 * x33;
                T x35 = x1 * x34;
                T x36 = x26 * x3;
                T x37 = x4 * z;
                T x38 = x36 * x37;
                T x39 = x1 * x38;
                T x40 = x11 * x12;
                T x41 = x36 * x40;
                T x42 = x15 * x18;
                T x43 = x26 * x29 * x42;
                T x44 = x1 * x26;
                T x45 = x42 * x44;
                T x46 = pow(x, 2);
                T x47 = x34 * x46;
                T x48 = x38 * x46;
                T x49 = x26 * x46;
                T x50 = x42 * x49;
                T x51 = pow(y, 2);
                T x52 = x27 * x51;
                T x53 = x29 * x52;
                T x54 = x37 * x51;
                T x55 = x49 * x54;
                T x56 = x26 * x40 * x51;
                T x57 = x44 * x54;
                T x58 = x26 * pow(z, 2);
                T x59 = x3 * x58;
                T x60 = x29 * x59;
                T x61 = x31 * x58;
                T x62 = x33 * x58;
                T x63 = x1 * x62;
                T x64 = x46 * x62;
                T x65 = x51 * x58;
                T x66 = x29 * x65;
                H[0] += x1 * x28;
                H[1] += x30;
                H[2] += x32;
                H[3] += x35;
                H[4] += x39;
                H[5] += x41;
                H[6] += x43;
                H[7] += x45;
                H[8] += x30;
                H[9] += x28 * x46;
                H[10] += x47;
                H[11] += x32;
                H[12] += x41;
                H[13] += x48;
                H[14] += x50;
                H[15] += x43;
                H[16] += x32;
                H[17] += x47;
                H[18] += x46 * x52;
                H[19] += x53;
                H[20] += x43;
                H[21] += x50;
                H[22] += x55;
                H[23] += x56;
                H[24] += x35;
                H[25] += x32;
                H[26] += x53;
                H[27] += x1 * x52;
                H[28] += x45;
                H[29] += x43;
                H[30] += x56;
                H[31] += x57;
                H[32] += x39;
                H[33] += x41;
                H[34] += x43;
                H[35] += x45;
                H[36] += x1 * x59;
                H[37] += x60;
                H[38] += x61;
                H[39] += x63;
                H[40] += x41;
                H[41] += x48;
                H[42] += x50;
                H[43] += x43;
                H[44] += x60;
                H[45] += x46 * x59;
                H[46] += x64;
                H[47] += x61;
                H[48] += x43;
                H[49] += x50;
                H[50] += x55;
                H[51] += x56;
                H[52] += x61;
                H[53] += x64;
                H[54] += x46 * x65;
                H[55] += x66;
                H[56] += x45;
                H[57] += x43;
                H[58] += x56;
                H[59] += x57;
                H[60] += x63;
                H[61] += x61;
                H[62] += x66;
                H[63] += x1 * x65;
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
                //	- Subexpressions: 36*ADD + 150*MUL + 42*SUB
                T x0 = 1.0 - x;
                T x1 = 1.0 - y;
                T x2 = 1.0 - z;
                T x3 = x1 * x2;
                T x4 = x0 * x3;
                T x5 = x * y;
                T x6 = x * x1;
                T x7 = x0 * y;
                T x8 = x0 * x1;
                T x9 = -px[0] * x8 - px[1] * x6 - px[2] * x5 - px[3] * x7 + px[4] * x0 * x1 + px[5] * x * x1 +
                       px[6] * x * y + px[7] * x0 * y;
                T x10 = x * z;
                T x11 = x * x2;
                T x12 = x0 * z;
                T x13 = x0 * x2;
                T x14 = -py[0] * x13 - py[1] * x11 + py[2] * x * x2 + py[3] * x0 * x2 - py[4] * x12 - py[5] * x10 +
                        py[6] * x * z + py[7] * x0 * z;
                T x15 = y * z;
                T x16 = x2 * y;
                T x17 = x1 * z;
                T x18 = -pz[0] * x3 + pz[1] * x1 * x2 + pz[2] * x2 * y - pz[3] * x16 - pz[4] * x17 + pz[5] * x1 * z +
                        pz[6] * y * z - pz[7] * x15;
                T x19 = -px[0] * x13 - px[1] * x11 + px[2] * x * x2 + px[3] * x0 * x2 - px[4] * x12 - px[5] * x10 +
                        px[6] * x * z + px[7] * x0 * z;
                T x20 = -py[0] * x3 + py[1] * x1 * x2 + py[2] * x2 * y - py[3] * x16 - py[4] * x17 + py[5] * x1 * z +
                        py[6] * y * z - py[7] * x15;
                T x21 = -pz[0] * x8 - pz[1] * x6 - pz[2] * x5 - pz[3] * x7 + pz[4] * x0 * x1 + pz[5] * x * x1 +
                        pz[6] * x * y + pz[7] * x0 * y;
                T x22 = -px[0] * x3 + px[1] * x1 * x2 + px[2] * x2 * y - px[3] * x16 - px[4] * x17 + px[5] * x1 * z +
                        px[6] * y * z - px[7] * x15;
                T x23 = -py[0] * x8 - py[1] * x6 - py[2] * x5 - py[3] * x7 + py[4] * x0 * x1 + py[5] * x * x1 +
                        py[6] * x * y + py[7] * x0 * y;
                T x24 = -pz[0] * x13 - pz[1] * x11 + pz[2] * x * x2 + pz[3] * x0 * x2 - pz[4] * x12 - pz[5] * x10 +
                        pz[6] * x * z + pz[7] * x0 * z;
                T x25 = x * x15;
                T x26 = x * x16;
                T x27 = x * x17;
                T x28 = x0 * x15;
                T x29 = x * x3;
                T x30 = x0 * x16;
                T x31 = x0 * x17;
                T x32 = weight *
                        (-x14 * x18 * x9 + x14 * x21 * x22 + x18 * x19 * x23 - x19 * x20 * x21 + x20 * x24 * x9 -
                         x22 * x23 * x24) *
                        (u[0] * x4 + u[1] * x29 + u[2] * x26 + u[3] * x30 + u[4] * x31 + u[5] * x27 + u[6] * x25 +
                         u[7] * x28);
                Hx[0] += x32 * x4;
                Hx[1] += x29 * x32;
                Hx[2] += x26 * x32;
                Hx[3] += x30 * x32;
                Hx[4] += x31 * x32;
                Hx[5] += x27 * x32;
                Hx[6] += x25 * x32;
                Hx[7] += x28 * x32;
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
                //	- Subexpressions: 36*ADD + 150*MUL + 42*SUB
                T x0 = 1.0 - x;
                T x1 = 1.0 - y;
                T x2 = 1.0 - z;
                T x3 = x1 * x2;
                T x4 = x0 * x3;
                T x5 = x * y;
                T x6 = x * x1;
                T x7 = x0 * y;
                T x8 = x0 * x1;
                T x9 = -px[0] * x8 - px[1] * x6 - px[2] * x5 - px[3] * x7 + px[4] * x0 * x1 + px[5] * x * x1 +
                       px[6] * x * y + px[7] * x0 * y;
                T x10 = x * z;
                T x11 = x * x2;
                T x12 = x0 * z;
                T x13 = x0 * x2;
                T x14 = -py[0] * x13 - py[1] * x11 + py[2] * x * x2 + py[3] * x0 * x2 - py[4] * x12 - py[5] * x10 +
                        py[6] * x * z + py[7] * x0 * z;
                T x15 = y * z;
                T x16 = x2 * y;
                T x17 = x1 * z;
                T x18 = -pz[0] * x3 + pz[1] * x1 * x2 + pz[2] * x2 * y - pz[3] * x16 - pz[4] * x17 + pz[5] * x1 * z +
                        pz[6] * y * z - pz[7] * x15;
                T x19 = -px[0] * x13 - px[1] * x11 + px[2] * x * x2 + px[3] * x0 * x2 - px[4] * x12 - px[5] * x10 +
                        px[6] * x * z + px[7] * x0 * z;
                T x20 = -py[0] * x3 + py[1] * x1 * x2 + py[2] * x2 * y - py[3] * x16 - py[4] * x17 + py[5] * x1 * z +
                        py[6] * y * z - py[7] * x15;
                T x21 = -pz[0] * x8 - pz[1] * x6 - pz[2] * x5 - pz[3] * x7 + pz[4] * x0 * x1 + pz[5] * x * x1 +
                        pz[6] * x * y + pz[7] * x0 * y;
                T x22 = -px[0] * x3 + px[1] * x1 * x2 + px[2] * x2 * y - px[3] * x16 - px[4] * x17 + px[5] * x1 * z +
                        px[6] * y * z - px[7] * x15;
                T x23 = -py[0] * x8 - py[1] * x6 - py[2] * x5 - py[3] * x7 + py[4] * x0 * x1 + py[5] * x * x1 +
                        py[6] * x * y + py[7] * x0 * y;
                T x24 = -pz[0] * x13 - pz[1] * x11 + pz[2] * x * x2 + pz[3] * x0 * x2 - pz[4] * x12 - pz[5] * x10 +
                        pz[6] * x * z + pz[7] * x0 * z;
                T x25 = x * x15;
                T x26 = x * x16;
                T x27 = x * x17;
                T x28 = x0 * x15;
                T x29 = x * x3;
                T x30 = x0 * x16;
                T x31 = x0 * x17;
                T x32 = weight *
                        (-x14 * x18 * x9 + x14 * x21 * x22 + x18 * x19 * x23 - x19 * x20 * x21 + x20 * x24 * x9 -
                         x22 * x23 * x24) *
                        (u[0] * x4 + u[1] * x29 + u[2] * x26 + u[3] * x30 + u[4] * x31 + u[5] * x27 + u[6] * x25 +
                         u[7] * x28);
                g[0] += x32 * x4;
                g[1] += x29 * x32;
                g[2] += x26 * x32;
                g[3] += x30 * x32;
                g[4] += x31 * x32;
                g[5] += x27 * x32;
                g[6] += x25 * x32;
                g[7] += x28 * x32;
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
                //	- Subexpressions: 36*ADD + 190*MUL + 6*POW + 42*SUB
                T x0 = 1.0 - x;
                T x1 = pow(x0, 2);
                T x2 = 1.0 - y;
                T x3 = pow(x2, 2);
                T x4 = 1.0 - z;
                T x5 = x * y;
                T x6 = x * x2;
                T x7 = x0 * y;
                T x8 = x0 * x2;
                T x9 = -px[0] * x8 - px[1] * x6 - px[2] * x5 - px[3] * x7 + px[4] * x0 * x2 + px[5] * x * x2 +
                       px[6] * x * y + px[7] * x0 * y;
                T x10 = x * z;
                T x11 = x * x4;
                T x12 = x0 * z;
                T x13 = x0 * x4;
                T x14 = -py[0] * x13 - py[1] * x11 + py[2] * x * x4 + py[3] * x0 * x4 - py[4] * x12 - py[5] * x10 +
                        py[6] * x * z + py[7] * x0 * z;
                T x15 = y * z;
                T x16 = x4 * y;
                T x17 = x2 * z;
                T x18 = x2 * x4;
                T x19 = -pz[0] * x18 + pz[1] * x2 * x4 + pz[2] * x4 * y - pz[3] * x16 - pz[4] * x17 + pz[5] * x2 * z +
                        pz[6] * y * z - pz[7] * x15;
                T x20 = -px[0] * x13 - px[1] * x11 + px[2] * x * x4 + px[3] * x0 * x4 - px[4] * x12 - px[5] * x10 +
                        px[6] * x * z + px[7] * x0 * z;
                T x21 = -py[0] * x18 + py[1] * x2 * x4 + py[2] * x4 * y - py[3] * x16 - py[4] * x17 + py[5] * x2 * z +
                        py[6] * y * z - py[7] * x15;
                T x22 = -pz[0] * x8 - pz[1] * x6 - pz[2] * x5 - pz[3] * x7 + pz[4] * x0 * x2 + pz[5] * x * x2 +
                        pz[6] * x * y + pz[7] * x0 * y;
                T x23 = -px[0] * x18 + px[1] * x2 * x4 + px[2] * x4 * y - px[3] * x16 - px[4] * x17 + px[5] * x2 * z +
                        px[6] * y * z - px[7] * x15;
                T x24 = -py[0] * x8 - py[1] * x6 - py[2] * x5 - py[3] * x7 + py[4] * x0 * x2 + py[5] * x * x2 +
                        py[6] * x * y + py[7] * x0 * y;
                T x25 = -pz[0] * x13 - pz[1] * x11 + pz[2] * x * x4 + pz[3] * x0 * x4 - pz[4] * x12 - pz[5] * x10 +
                        pz[6] * x * z + pz[7] * x0 * z;
                T x26 = weight * (-x14 * x19 * x9 + x14 * x22 * x23 + x19 * x20 * x24 - x20 * x21 * x22 +
                                  x21 * x25 * x9 - x23 * x24 * x25);
                T x27 = x26 * pow(x4, 2);
                T x28 = x27 * x3;
                T x29 = x * x0;
                T x30 = x28 * x29;
                T x31 = x5 * x8;
                T x32 = x27 * x31;
                T x33 = x2 * y;
                T x34 = x27 * x33;
                T x35 = x1 * x34;
                T x36 = x26 * x3;
                T x37 = x4 * z;
                T x38 = x36 * x37;
                T x39 = x1 * x38;
                T x40 = x11 * x12;
                T x41 = x36 * x40;
                T x42 = x0 * x18;
                T x43 = x * x15;
                T x44 = x26 * x42 * x43;
                T x45 = x1 * x26;
                T x46 = x15 * x18;
                T x47 = x45 * x46;
                T x48 = pow(x, 2);
                T x49 = x34 * x48;
                T x50 = x38 * x48;
                T x51 = x26 * x48;
                T x52 = x46 * x51;
                T x53 = pow(y, 2);
                T x54 = x27 * x53;
                T x55 = x29 * x54;
                T x56 = x37 * x53;
                T x57 = x51 * x56;
                T x58 = x26 * x40 * x53;
                T x59 = x45 * x56;
                T x60 = x26 * pow(z, 2);
                T x61 = x3 * x60;
                T x62 = x29 * x61;
                T x63 = x31 * x60;
                T x64 = x33 * x60;
                T x65 = x1 * x64;
                T x66 = x48 * x64;
                T x67 = x53 * x60;
                T x68 = x29 * x67;
                T x69 = x * x16;
                T x70 = x * x17;
                T x71 = x0 * x15;
                T x72 = x * x18;
                T x73 = x0 * x16;
                T x74 = x0 * x17;
                T x75 = x26 * (u[0] * x42 + u[1] * x72 + u[2] * x69 + u[3] * x73 + u[4] * x74 + u[5] * x70 +
                               u[6] * x43 + u[7] * x71);
                H[0] += x1 * x28;
                H[1] += x30;
                H[2] += x32;
                H[3] += x35;
                H[4] += x39;
                H[5] += x41;
                H[6] += x44;
                H[7] += x47;
                H[8] += x30;
                H[9] += x28 * x48;
                H[10] += x49;
                H[11] += x32;
                H[12] += x41;
                H[13] += x50;
                H[14] += x52;
                H[15] += x44;
                H[16] += x32;
                H[17] += x49;
                H[18] += x48 * x54;
                H[19] += x55;
                H[20] += x44;
                H[21] += x52;
                H[22] += x57;
                H[23] += x58;
                H[24] += x35;
                H[25] += x32;
                H[26] += x55;
                H[27] += x1 * x54;
                H[28] += x47;
                H[29] += x44;
                H[30] += x58;
                H[31] += x59;
                H[32] += x39;
                H[33] += x41;
                H[34] += x44;
                H[35] += x47;
                H[36] += x1 * x61;
                H[37] += x62;
                H[38] += x63;
                H[39] += x65;
                H[40] += x41;
                H[41] += x50;
                H[42] += x52;
                H[43] += x44;
                H[44] += x62;
                H[45] += x48 * x61;
                H[46] += x66;
                H[47] += x63;
                H[48] += x44;
                H[49] += x52;
                H[50] += x57;
                H[51] += x58;
                H[52] += x63;
                H[53] += x66;
                H[54] += x48 * x67;
                H[55] += x68;
                H[56] += x47;
                H[57] += x44;
                H[58] += x58;
                H[59] += x59;
                H[60] += x65;
                H[61] += x63;
                H[62] += x68;
                H[63] += x1 * x67;
                g[0] += x42 * x75;
                g[1] += x72 * x75;
                g[2] += x69 * x75;
                g[3] += x73 * x75;
                g[4] += x74 * x75;
                g[5] += x70 * x75;
                g[6] += x43 * x75;
                g[7] += x71 * x75;
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using MassHex8 = utopia::kokkos::
            AutoKernel<FE, utopia::kernels::Mass<utopia::kernels::Hex8<typename FE::Scalar, typename FE::Scalar>>, 3>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_Mass_3_IMPL_hpp
