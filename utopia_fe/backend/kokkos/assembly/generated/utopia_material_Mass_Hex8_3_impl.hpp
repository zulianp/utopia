#ifndef UTOPIA_TPL_MATERIAL_Mass_Hex8_3_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_Mass_Hex8_3_IMPL_hpp

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
                //	- Subexpressions: 35*ADD + 46*DIV + 174*MUL + 6*POW + 45*SUB
                T x0 = pow(1 - z, 2);
                T x1 = (1.0 / 2.0) * x;
                T x2 = 0.5 - x1;
                T x3 = (1.0 / 2.0) * y;
                T x4 = 0.5 - x3;
                T x5 = (1.0 / 2.0) * x4;
                T x6 = x2 * x5;
                T x7 = x1 + 0.5;
                T x8 = x5 * x7;
                T x9 = x3 + 0.5;
                T x10 = (1.0 / 2.0) * x9;
                T x11 = x10 * x7;
                T x12 = x10 * x2;
                T x13 = -px[0] * x6 - px[1] * x8 - px[2] * x11 - px[3] * x12 + (1.0 / 2.0) * px[4] * x2 * x4 +
                        (1.0 / 2.0) * px[5] * x4 * x7 + (1.0 / 2.0) * px[6] * x7 * x9 + (1.0 / 2.0) * px[7] * x2 * x9;
                T x14 = (1.0 / 2.0) * z;
                T x15 = 0.5 - x14;
                T x16 = (1.0 / 2.0) * x15;
                T x17 = x16 * x2;
                T x18 = x16 * x7;
                T x19 = x14 + 0.5;
                T x20 = (1.0 / 2.0) * x19;
                T x21 = x2 * x20;
                T x22 = x20 * x7;
                T x23 = -py[0] * x17 - py[1] * x18 + (1.0 / 2.0) * py[2] * x15 * x7 + (1.0 / 2.0) * py[3] * x15 * x2 -
                        py[4] * x21 - py[5] * x22 + (1.0 / 2.0) * py[6] * x19 * x7 + (1.0 / 2.0) * py[7] * x19 * x2;
                T x24 = x16 * x4;
                T x25 = x16 * x9;
                T x26 = x20 * x4;
                T x27 = x20 * x9;
                T x28 = -pz[0] * x24 + (1.0 / 2.0) * pz[1] * x15 * x4 + (1.0 / 2.0) * pz[2] * x15 * x9 - pz[3] * x25 -
                        pz[4] * x26 + (1.0 / 2.0) * pz[5] * x19 * x4 + (1.0 / 2.0) * pz[6] * x19 * x9 - pz[7] * x27;
                T x29 = -px[0] * x17 - px[1] * x18 + (1.0 / 2.0) * px[2] * x15 * x7 + (1.0 / 2.0) * px[3] * x15 * x2 -
                        px[4] * x21 - px[5] * x22 + (1.0 / 2.0) * px[6] * x19 * x7 + (1.0 / 2.0) * px[7] * x19 * x2;
                T x30 = -py[0] * x24 + (1.0 / 2.0) * py[1] * x15 * x4 + (1.0 / 2.0) * py[2] * x15 * x9 - py[3] * x25 -
                        py[4] * x26 + (1.0 / 2.0) * py[5] * x19 * x4 + (1.0 / 2.0) * py[6] * x19 * x9 - py[7] * x27;
                T x31 = -pz[0] * x6 - pz[1] * x8 - pz[2] * x11 - pz[3] * x12 + (1.0 / 2.0) * pz[4] * x2 * x4 +
                        (1.0 / 2.0) * pz[5] * x4 * x7 + (1.0 / 2.0) * pz[6] * x7 * x9 + (1.0 / 2.0) * pz[7] * x2 * x9;
                T x32 = -px[0] * x24 + (1.0 / 2.0) * px[1] * x15 * x4 + (1.0 / 2.0) * px[2] * x15 * x9 - px[3] * x25 -
                        px[4] * x26 + (1.0 / 2.0) * px[5] * x19 * x4 + (1.0 / 2.0) * px[6] * x19 * x9 - px[7] * x27;
                T x33 = -py[0] * x6 - py[1] * x8 - py[2] * x11 - py[3] * x12 + (1.0 / 2.0) * py[4] * x2 * x4 +
                        (1.0 / 2.0) * py[5] * x4 * x7 + (1.0 / 2.0) * py[6] * x7 * x9 + (1.0 / 2.0) * py[7] * x2 * x9;
                T x34 = -pz[0] * x17 - pz[1] * x18 + (1.0 / 2.0) * pz[2] * x15 * x7 + (1.0 / 2.0) * pz[3] * x15 * x2 -
                        pz[4] * x21 - pz[5] * x22 + (1.0 / 2.0) * pz[6] * x19 * x7 + (1.0 / 2.0) * pz[7] * x19 * x2;
                T x35 = weight * (-x13 * x23 * x28 + x13 * x30 * x34 + x23 * x31 * x32 + x28 * x29 * x33 -
                                  x29 * x30 * x31 - x32 * x33 * x34);
                T x36 = (1.0 / 64.0) * x35;
                T x37 = x0 * x36;
                T x38 = pow(1 - x, 2);
                T x39 = pow(1 - y, 2);
                T x40 = x38 * x39;
                T x41 = (1.0 / 16.0) * x35;
                T x42 = x0 * x41;
                T x43 = x2 * x7;
                T x44 = x39 * x43;
                T x45 = x42 * x44;
                T x46 = (1.0 / 4.0) * x35;
                T x47 = x4 * x9;
                T x48 = x43 * x47;
                T x49 = x46 * x48;
                T x50 = x0 * x49;
                T x51 = x42 * x47;
                T x52 = x38 * x51;
                T x53 = x15 * x19;
                T x54 = x41 * x53;
                T x55 = x40 * x54;
                T x56 = x46 * x53;
                T x57 = x44 * x56;
                T x58 = x35 * x48 * x53;
                T x59 = x47 * x56;
                T x60 = x38 * x59;
                T x61 = pow(x + 1, 2);
                T x62 = x37 * x61;
                T x63 = x51 * x61;
                T x64 = x54 * x61;
                T x65 = x39 * x64;
                T x66 = x59 * x61;
                T x67 = pow(y + 1, 2);
                T x68 = x43 * x67;
                T x69 = x42 * x68;
                T x70 = x64 * x67;
                T x71 = x56 * x68;
                T x72 = x38 * x67;
                T x73 = x54 * x72;
                T x74 = pow(z + 1, 2);
                T x75 = x36 * x74;
                T x76 = x41 * x74;
                T x77 = x44 * x76;
                T x78 = x49 * x74;
                T x79 = x47 * x76;
                T x80 = x38 * x79;
                T x81 = x61 * x75;
                T x82 = x61 * x79;
                T x83 = x68 * x76;
                H[0] += x37 * x40;
                H[1] += x45;
                H[2] += x50;
                H[3] += x52;
                H[4] += x55;
                H[5] += x57;
                H[6] += x58;
                H[7] += x60;
                H[8] += x45;
                H[9] += x39 * x62;
                H[10] += x63;
                H[11] += x50;
                H[12] += x57;
                H[13] += x65;
                H[14] += x66;
                H[15] += x58;
                H[16] += x50;
                H[17] += x63;
                H[18] += x62 * x67;
                H[19] += x69;
                H[20] += x58;
                H[21] += x66;
                H[22] += x70;
                H[23] += x71;
                H[24] += x52;
                H[25] += x50;
                H[26] += x69;
                H[27] += x37 * x72;
                H[28] += x60;
                H[29] += x58;
                H[30] += x71;
                H[31] += x73;
                H[32] += x55;
                H[33] += x57;
                H[34] += x58;
                H[35] += x60;
                H[36] += x40 * x75;
                H[37] += x77;
                H[38] += x78;
                H[39] += x80;
                H[40] += x57;
                H[41] += x65;
                H[42] += x66;
                H[43] += x58;
                H[44] += x77;
                H[45] += x39 * x81;
                H[46] += x82;
                H[47] += x78;
                H[48] += x58;
                H[49] += x66;
                H[50] += x70;
                H[51] += x71;
                H[52] += x78;
                H[53] += x82;
                H[54] += x67 * x81;
                H[55] += x83;
                H[56] += x60;
                H[57] += x58;
                H[58] += x71;
                H[59] += x73;
                H[60] += x80;
                H[61] += x78;
                H[62] += x83;
                H[63] += x72 * x75;
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
                //	- Subexpressions: 39*ADD + 43*DIV + 154*MUL + 42*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * z;
                T x5 = 0.5 - x4;
                T x6 = x3 * x5;
                T x7 = x1 * x6;
                T x8 = (1.0 / 2.0) * x3;
                T x9 = x1 * x8;
                T x10 = x0 + 0.5;
                T x11 = x10 * x8;
                T x12 = x2 + 0.5;
                T x13 = (1.0 / 2.0) * x12;
                T x14 = x10 * x13;
                T x15 = x1 * x13;
                T x16 = -px[0] * x9 - px[1] * x11 - px[2] * x14 - px[3] * x15 + (1.0 / 2.0) * px[4] * x1 * x3 +
                        (1.0 / 2.0) * px[5] * x10 * x3 + (1.0 / 2.0) * px[6] * x10 * x12 +
                        (1.0 / 2.0) * px[7] * x1 * x12;
                T x17 = (1.0 / 2.0) * x5;
                T x18 = x1 * x17;
                T x19 = x10 * x17;
                T x20 = x4 + 0.5;
                T x21 = (1.0 / 2.0) * x20;
                T x22 = x1 * x21;
                T x23 = x10 * x21;
                T x24 = -py[0] * x18 - py[1] * x19 + (1.0 / 2.0) * py[2] * x10 * x5 + (1.0 / 2.0) * py[3] * x1 * x5 -
                        py[4] * x22 - py[5] * x23 + (1.0 / 2.0) * py[6] * x10 * x20 + (1.0 / 2.0) * py[7] * x1 * x20;
                T x25 = x17 * x3;
                T x26 = x12 * x17;
                T x27 = x21 * x3;
                T x28 = x12 * x21;
                T x29 = -pz[0] * x25 + (1.0 / 2.0) * pz[1] * x3 * x5 + (1.0 / 2.0) * pz[2] * x12 * x5 - pz[3] * x26 -
                        pz[4] * x27 + (1.0 / 2.0) * pz[5] * x20 * x3 + (1.0 / 2.0) * pz[6] * x12 * x20 - pz[7] * x28;
                T x30 = -px[0] * x18 - px[1] * x19 + (1.0 / 2.0) * px[2] * x10 * x5 + (1.0 / 2.0) * px[3] * x1 * x5 -
                        px[4] * x22 - px[5] * x23 + (1.0 / 2.0) * px[6] * x10 * x20 + (1.0 / 2.0) * px[7] * x1 * x20;
                T x31 = -py[0] * x25 + (1.0 / 2.0) * py[1] * x3 * x5 + (1.0 / 2.0) * py[2] * x12 * x5 - py[3] * x26 -
                        py[4] * x27 + (1.0 / 2.0) * py[5] * x20 * x3 + (1.0 / 2.0) * py[6] * x12 * x20 - py[7] * x28;
                T x32 = -pz[0] * x9 - pz[1] * x11 - pz[2] * x14 - pz[3] * x15 + (1.0 / 2.0) * pz[4] * x1 * x3 +
                        (1.0 / 2.0) * pz[5] * x10 * x3 + (1.0 / 2.0) * pz[6] * x10 * x12 +
                        (1.0 / 2.0) * pz[7] * x1 * x12;
                T x33 = -px[0] * x25 + (1.0 / 2.0) * px[1] * x3 * x5 + (1.0 / 2.0) * px[2] * x12 * x5 - px[3] * x26 -
                        px[4] * x27 + (1.0 / 2.0) * px[5] * x20 * x3 + (1.0 / 2.0) * px[6] * x12 * x20 - px[7] * x28;
                T x34 = -py[0] * x9 - py[1] * x11 - py[2] * x14 - py[3] * x15 + (1.0 / 2.0) * py[4] * x1 * x3 +
                        (1.0 / 2.0) * py[5] * x10 * x3 + (1.0 / 2.0) * py[6] * x10 * x12 +
                        (1.0 / 2.0) * py[7] * x1 * x12;
                T x35 = -pz[0] * x18 - pz[1] * x19 + (1.0 / 2.0) * pz[2] * x10 * x5 + (1.0 / 2.0) * pz[3] * x1 * x5 -
                        pz[4] * x22 - pz[5] * x23 + (1.0 / 2.0) * pz[6] * x10 * x20 + (1.0 / 2.0) * pz[7] * x1 * x20;
                T x36 = x10 * x6;
                T x37 = x12 * x5;
                T x38 = x10 * x37;
                T x39 = x1 * x37;
                T x40 = x20 * x3;
                T x41 = x1 * x40;
                T x42 = x10 * x40;
                T x43 = x12 * x20;
                T x44 = x10 * x43;
                T x45 = x1 * x43;
                T x46 = weight *
                        (-x16 * x24 * x29 + x16 * x31 * x35 + x24 * x32 * x33 + x29 * x30 * x34 - x30 * x31 * x32 -
                         x33 * x34 * x35) *
                        (u[0] * x7 + u[1] * x36 + u[2] * x38 + u[3] * x39 + u[4] * x41 + u[5] * x42 + u[6] * x44 +
                         u[7] * x45);
                Hx[0] += x46 * x7;
                Hx[1] += x36 * x46;
                Hx[2] += x38 * x46;
                Hx[3] += x39 * x46;
                Hx[4] += x41 * x46;
                Hx[5] += x42 * x46;
                Hx[6] += x44 * x46;
                Hx[7] += x45 * x46;
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
                //	- Subexpressions: 39*ADD + 43*DIV + 154*MUL + 42*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * z;
                T x5 = 0.5 - x4;
                T x6 = x3 * x5;
                T x7 = x1 * x6;
                T x8 = (1.0 / 2.0) * x3;
                T x9 = x1 * x8;
                T x10 = x0 + 0.5;
                T x11 = x10 * x8;
                T x12 = x2 + 0.5;
                T x13 = (1.0 / 2.0) * x12;
                T x14 = x10 * x13;
                T x15 = x1 * x13;
                T x16 = -px[0] * x9 - px[1] * x11 - px[2] * x14 - px[3] * x15 + (1.0 / 2.0) * px[4] * x1 * x3 +
                        (1.0 / 2.0) * px[5] * x10 * x3 + (1.0 / 2.0) * px[6] * x10 * x12 +
                        (1.0 / 2.0) * px[7] * x1 * x12;
                T x17 = (1.0 / 2.0) * x5;
                T x18 = x1 * x17;
                T x19 = x10 * x17;
                T x20 = x4 + 0.5;
                T x21 = (1.0 / 2.0) * x20;
                T x22 = x1 * x21;
                T x23 = x10 * x21;
                T x24 = -py[0] * x18 - py[1] * x19 + (1.0 / 2.0) * py[2] * x10 * x5 + (1.0 / 2.0) * py[3] * x1 * x5 -
                        py[4] * x22 - py[5] * x23 + (1.0 / 2.0) * py[6] * x10 * x20 + (1.0 / 2.0) * py[7] * x1 * x20;
                T x25 = x17 * x3;
                T x26 = x12 * x17;
                T x27 = x21 * x3;
                T x28 = x12 * x21;
                T x29 = -pz[0] * x25 + (1.0 / 2.0) * pz[1] * x3 * x5 + (1.0 / 2.0) * pz[2] * x12 * x5 - pz[3] * x26 -
                        pz[4] * x27 + (1.0 / 2.0) * pz[5] * x20 * x3 + (1.0 / 2.0) * pz[6] * x12 * x20 - pz[7] * x28;
                T x30 = -px[0] * x18 - px[1] * x19 + (1.0 / 2.0) * px[2] * x10 * x5 + (1.0 / 2.0) * px[3] * x1 * x5 -
                        px[4] * x22 - px[5] * x23 + (1.0 / 2.0) * px[6] * x10 * x20 + (1.0 / 2.0) * px[7] * x1 * x20;
                T x31 = -py[0] * x25 + (1.0 / 2.0) * py[1] * x3 * x5 + (1.0 / 2.0) * py[2] * x12 * x5 - py[3] * x26 -
                        py[4] * x27 + (1.0 / 2.0) * py[5] * x20 * x3 + (1.0 / 2.0) * py[6] * x12 * x20 - py[7] * x28;
                T x32 = -pz[0] * x9 - pz[1] * x11 - pz[2] * x14 - pz[3] * x15 + (1.0 / 2.0) * pz[4] * x1 * x3 +
                        (1.0 / 2.0) * pz[5] * x10 * x3 + (1.0 / 2.0) * pz[6] * x10 * x12 +
                        (1.0 / 2.0) * pz[7] * x1 * x12;
                T x33 = -px[0] * x25 + (1.0 / 2.0) * px[1] * x3 * x5 + (1.0 / 2.0) * px[2] * x12 * x5 - px[3] * x26 -
                        px[4] * x27 + (1.0 / 2.0) * px[5] * x20 * x3 + (1.0 / 2.0) * px[6] * x12 * x20 - px[7] * x28;
                T x34 = -py[0] * x9 - py[1] * x11 - py[2] * x14 - py[3] * x15 + (1.0 / 2.0) * py[4] * x1 * x3 +
                        (1.0 / 2.0) * py[5] * x10 * x3 + (1.0 / 2.0) * py[6] * x10 * x12 +
                        (1.0 / 2.0) * py[7] * x1 * x12;
                T x35 = -pz[0] * x18 - pz[1] * x19 + (1.0 / 2.0) * pz[2] * x10 * x5 + (1.0 / 2.0) * pz[3] * x1 * x5 -
                        pz[4] * x22 - pz[5] * x23 + (1.0 / 2.0) * pz[6] * x10 * x20 + (1.0 / 2.0) * pz[7] * x1 * x20;
                T x36 = x10 * x6;
                T x37 = x12 * x5;
                T x38 = x10 * x37;
                T x39 = x1 * x37;
                T x40 = x20 * x3;
                T x41 = x1 * x40;
                T x42 = x10 * x40;
                T x43 = x12 * x20;
                T x44 = x10 * x43;
                T x45 = x1 * x43;
                T x46 = weight *
                        (-x16 * x24 * x29 + x16 * x31 * x35 + x24 * x32 * x33 + x29 * x30 * x34 - x30 * x31 * x32 -
                         x33 * x34 * x35) *
                        (u[0] * x7 + u[1] * x36 + u[2] * x38 + u[3] * x39 + u[4] * x41 + u[5] * x42 + u[6] * x44 +
                         u[7] * x45);
                g[0] += x46 * x7;
                g[1] += x36 * x46;
                g[2] += x38 * x46;
                g[3] += x39 * x46;
                g[4] += x41 * x46;
                g[5] += x42 * x46;
                g[6] += x44 * x46;
                g[7] += x45 * x46;
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
                //	- Subexpressions: 42*ADD + 46*DIV + 196*MUL + 6*POW + 45*SUB
                T x0 = pow(1 - z, 2);
                T x1 = (1.0 / 2.0) * x;
                T x2 = 0.5 - x1;
                T x3 = (1.0 / 2.0) * y;
                T x4 = 0.5 - x3;
                T x5 = (1.0 / 2.0) * x4;
                T x6 = x2 * x5;
                T x7 = x1 + 0.5;
                T x8 = x5 * x7;
                T x9 = x3 + 0.5;
                T x10 = (1.0 / 2.0) * x9;
                T x11 = x10 * x7;
                T x12 = x10 * x2;
                T x13 = -px[0] * x6 - px[1] * x8 - px[2] * x11 - px[3] * x12 + (1.0 / 2.0) * px[4] * x2 * x4 +
                        (1.0 / 2.0) * px[5] * x4 * x7 + (1.0 / 2.0) * px[6] * x7 * x9 + (1.0 / 2.0) * px[7] * x2 * x9;
                T x14 = (1.0 / 2.0) * z;
                T x15 = 0.5 - x14;
                T x16 = (1.0 / 2.0) * x15;
                T x17 = x16 * x2;
                T x18 = x16 * x7;
                T x19 = x14 + 0.5;
                T x20 = (1.0 / 2.0) * x19;
                T x21 = x2 * x20;
                T x22 = x20 * x7;
                T x23 = -py[0] * x17 - py[1] * x18 + (1.0 / 2.0) * py[2] * x15 * x7 + (1.0 / 2.0) * py[3] * x15 * x2 -
                        py[4] * x21 - py[5] * x22 + (1.0 / 2.0) * py[6] * x19 * x7 + (1.0 / 2.0) * py[7] * x19 * x2;
                T x24 = x16 * x4;
                T x25 = x16 * x9;
                T x26 = x20 * x4;
                T x27 = x20 * x9;
                T x28 = -pz[0] * x24 + (1.0 / 2.0) * pz[1] * x15 * x4 + (1.0 / 2.0) * pz[2] * x15 * x9 - pz[3] * x25 -
                        pz[4] * x26 + (1.0 / 2.0) * pz[5] * x19 * x4 + (1.0 / 2.0) * pz[6] * x19 * x9 - pz[7] * x27;
                T x29 = -px[0] * x17 - px[1] * x18 + (1.0 / 2.0) * px[2] * x15 * x7 + (1.0 / 2.0) * px[3] * x15 * x2 -
                        px[4] * x21 - px[5] * x22 + (1.0 / 2.0) * px[6] * x19 * x7 + (1.0 / 2.0) * px[7] * x19 * x2;
                T x30 = -py[0] * x24 + (1.0 / 2.0) * py[1] * x15 * x4 + (1.0 / 2.0) * py[2] * x15 * x9 - py[3] * x25 -
                        py[4] * x26 + (1.0 / 2.0) * py[5] * x19 * x4 + (1.0 / 2.0) * py[6] * x19 * x9 - py[7] * x27;
                T x31 = -pz[0] * x6 - pz[1] * x8 - pz[2] * x11 - pz[3] * x12 + (1.0 / 2.0) * pz[4] * x2 * x4 +
                        (1.0 / 2.0) * pz[5] * x4 * x7 + (1.0 / 2.0) * pz[6] * x7 * x9 + (1.0 / 2.0) * pz[7] * x2 * x9;
                T x32 = -px[0] * x24 + (1.0 / 2.0) * px[1] * x15 * x4 + (1.0 / 2.0) * px[2] * x15 * x9 - px[3] * x25 -
                        px[4] * x26 + (1.0 / 2.0) * px[5] * x19 * x4 + (1.0 / 2.0) * px[6] * x19 * x9 - px[7] * x27;
                T x33 = -py[0] * x6 - py[1] * x8 - py[2] * x11 - py[3] * x12 + (1.0 / 2.0) * py[4] * x2 * x4 +
                        (1.0 / 2.0) * py[5] * x4 * x7 + (1.0 / 2.0) * py[6] * x7 * x9 + (1.0 / 2.0) * py[7] * x2 * x9;
                T x34 = -pz[0] * x17 - pz[1] * x18 + (1.0 / 2.0) * pz[2] * x15 * x7 + (1.0 / 2.0) * pz[3] * x15 * x2 -
                        pz[4] * x21 - pz[5] * x22 + (1.0 / 2.0) * pz[6] * x19 * x7 + (1.0 / 2.0) * pz[7] * x19 * x2;
                T x35 = weight * (-x13 * x23 * x28 + x13 * x30 * x34 + x23 * x31 * x32 + x28 * x29 * x33 -
                                  x29 * x30 * x31 - x32 * x33 * x34);
                T x36 = (1.0 / 64.0) * x35;
                T x37 = x0 * x36;
                T x38 = pow(1 - x, 2);
                T x39 = pow(1 - y, 2);
                T x40 = x38 * x39;
                T x41 = (1.0 / 16.0) * x35;
                T x42 = x0 * x41;
                T x43 = x2 * x7;
                T x44 = x39 * x43;
                T x45 = x42 * x44;
                T x46 = x4 * x9;
                T x47 = (1.0 / 4.0) * x35;
                T x48 = x43 * x46 * x47;
                T x49 = x0 * x48;
                T x50 = x42 * x46;
                T x51 = x38 * x50;
                T x52 = x15 * x19;
                T x53 = x41 * x52;
                T x54 = x40 * x53;
                T x55 = x47 * x52;
                T x56 = x44 * x55;
                T x57 = x15 * x4;
                T x58 = x2 * x57;
                T x59 = x19 * x9;
                T x60 = x59 * x7;
                T x61 = x35 * x58 * x60;
                T x62 = x47 * x57 * x59;
                T x63 = x38 * x62;
                T x64 = pow(x + 1, 2);
                T x65 = x37 * x64;
                T x66 = x50 * x64;
                T x67 = x53 * x64;
                T x68 = x39 * x67;
                T x69 = x62 * x64;
                T x70 = pow(y + 1, 2);
                T x71 = x43 * x70;
                T x72 = x42 * x71;
                T x73 = x67 * x70;
                T x74 = x55 * x71;
                T x75 = x38 * x70;
                T x76 = x53 * x75;
                T x77 = pow(z + 1, 2);
                T x78 = x36 * x77;
                T x79 = x41 * x77;
                T x80 = x44 * x79;
                T x81 = x48 * x77;
                T x82 = x46 * x79;
                T x83 = x38 * x82;
                T x84 = x64 * x78;
                T x85 = x64 * x82;
                T x86 = x71 * x79;
                T x87 = x57 * x7;
                T x88 = x15 * x9;
                T x89 = x7 * x88;
                T x90 = x2 * x88;
                T x91 = x19 * x4;
                T x92 = x2 * x91;
                T x93 = x7 * x91;
                T x94 = x2 * x59;
                T x95 = x35 * (u[0] * x58 + u[1] * x87 + u[2] * x89 + u[3] * x90 + u[4] * x92 + u[5] * x93 +
                               u[6] * x60 + u[7] * x94);
                H[0] += x37 * x40;
                H[1] += x45;
                H[2] += x49;
                H[3] += x51;
                H[4] += x54;
                H[5] += x56;
                H[6] += x61;
                H[7] += x63;
                H[8] += x45;
                H[9] += x39 * x65;
                H[10] += x66;
                H[11] += x49;
                H[12] += x56;
                H[13] += x68;
                H[14] += x69;
                H[15] += x61;
                H[16] += x49;
                H[17] += x66;
                H[18] += x65 * x70;
                H[19] += x72;
                H[20] += x61;
                H[21] += x69;
                H[22] += x73;
                H[23] += x74;
                H[24] += x51;
                H[25] += x49;
                H[26] += x72;
                H[27] += x37 * x75;
                H[28] += x63;
                H[29] += x61;
                H[30] += x74;
                H[31] += x76;
                H[32] += x54;
                H[33] += x56;
                H[34] += x61;
                H[35] += x63;
                H[36] += x40 * x78;
                H[37] += x80;
                H[38] += x81;
                H[39] += x83;
                H[40] += x56;
                H[41] += x68;
                H[42] += x69;
                H[43] += x61;
                H[44] += x80;
                H[45] += x39 * x84;
                H[46] += x85;
                H[47] += x81;
                H[48] += x61;
                H[49] += x69;
                H[50] += x73;
                H[51] += x74;
                H[52] += x81;
                H[53] += x85;
                H[54] += x70 * x84;
                H[55] += x86;
                H[56] += x63;
                H[57] += x61;
                H[58] += x74;
                H[59] += x76;
                H[60] += x83;
                H[61] += x81;
                H[62] += x86;
                H[63] += x75 * x78;
                g[0] += x58 * x95;
                g[1] += x87 * x95;
                g[2] += x89 * x95;
                g[3] += x90 * x95;
                g[4] += x92 * x95;
                g[5] += x93 * x95;
                g[6] += x60 * x95;
                g[7] += x94 * x95;
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

#endif  // UTOPIA_TPL_MATERIAL_Mass_Hex8_3_IMPL_hpp
