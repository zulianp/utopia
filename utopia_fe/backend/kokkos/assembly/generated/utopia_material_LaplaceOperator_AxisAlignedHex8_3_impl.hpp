#ifndef UTOPIA_TPL_MATERIAL_LaplaceOperator_AxisAlignedHex8_3_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_LaplaceOperator_AxisAlignedHex8_3_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_AxisAlignedHex8_3.hpp"
#include "utopia_material_LaplaceOperator.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of LaplaceOperator for symmetric element pair trial=test=AxisAlignedHex8
         */
        template <typename T, typename GeoT>
        class LaplaceOperator<AxisAlignedHex8<T, GeoT>> {
        public:
            using ElemT = AxisAlignedHex8<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "LaplaceOperator<AxisAlignedHex8>"; }

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
                //	- Result: 8*ADD + 64*ADDAUGMENTEDASSIGNMENT + 32*MUL
                //	- Subexpressions: 27*ADD + 10*DIV + 127*MUL + 19*NEG + 12*POW + 35*SUB
                T x0 = px[0] - px[6];
                T x1 = (1.0 / 64.0) / pow(x0, 2);
                T x2 = y - 1.0;
                T x3 = pow(x2, 2);
                T x4 = z - 1.0;
                T x5 = pow(x4, 2);
                T x6 = x3 * x5;
                T x7 = x - 1.0;
                T x8 = pow(x7, 2);
                T x9 = (1.0 / 64.0) * x8;
                T x10 = py[0] - py[6];
                T x11 = pow(x10, -2);
                T x12 = x11 * x5;
                T x13 = pz[0] - pz[6];
                T x14 = pow(x13, -2);
                T x15 = x14 * x3;
                T x16 = -weight * x0 * x10 * x13;
                T x17 = 8 * px[0] - 8 * px[6];
                T x18 = (1.0 / 8.0) / (x0 * x17);
                T x19 = 8 * py[0] - 8 * py[6];
                T x20 = 1 / (x10 * x19);
                T x21 = x20 * x5;
                T x22 = x + 1.0;
                T x23 = x22 * x7;
                T x24 = (1.0 / 8.0) * x23;
                T x25 = x21 * x24;
                T x26 = 8 * pz[0] - 8 * pz[6];
                T x27 = 1 / (x13 * x26);
                T x28 = x27 * x3;
                T x29 = x24 * x28;
                T x30 = x16 * (-x18 * x6 - x25 - x29);
                T x31 = y + 1.0;
                T x32 = x2 * x31;
                T x33 = x1 * x32;
                T x34 = (1.0 / 64.0) * x23;
                T x35 = x14 * x32;
                T x36 = x34 * x35;
                T x37 = x16 * (x12 * x34 + x33 * x5 + x36);
                T x38 = (1.0 / 8.0) * x8;
                T x39 = x18 * x32;
                T x40 = x39 * x5;
                T x41 = x27 * x32;
                T x42 = x38 * x41;
                T x43 = x16 * (-x21 * x38 - x40 - x42);
                T x44 = z + 1.0;
                T x45 = x4 * x44;
                T x46 = x18 * x3;
                T x47 = x45 * x46;
                T x48 = x20 * x45;
                T x49 = x38 * x48;
                T x50 = x16 * (-x28 * x38 - x47 - x49);
                T x51 = x1 * x3;
                T x52 = x11 * x45;
                T x53 = x34 * x52;
                T x54 = x16 * (x15 * x34 + x45 * x51 + x53);
                T x55 = x16 * (-x24 * x41 - x24 * x48 - x39 * x45);
                T x56 = x33 * x45;
                T x57 = x16 * (x35 * x9 + x52 * x9 + x56);
                T x58 = pow(x17, -2);
                T x59 = x5 * x58;
                T x60 = pow(x22, 2);
                T x61 = pow(x19, -2);
                T x62 = x5 * x61;
                T x63 = pow(x26, -2);
                T x64 = x3 * x63;
                T x65 = (1.0 / 8.0) * x60;
                T x66 = x41 * x65;
                T x67 = x16 * (-x21 * x65 - x40 - x66);
                T x68 = x32 * x63;
                T x69 = x23 * x68;
                T x70 = x16 * (x23 * x62 + x32 * x59 + x69);
                T x71 = x45 * x58;
                T x72 = x45 * x61;
                T x73 = x23 * x72;
                T x74 = x16 * (x23 * x64 + x3 * x71 + x73);
                T x75 = x48 * x65;
                T x76 = x16 * (-x28 * x65 - x47 - x75);
                T x77 = x32 * x71;
                T x78 = x16 * (x60 * x68 + x60 * x72 + x77);
                T x79 = pow(x31, 2);
                T x80 = x1 * x79;
                T x81 = (1.0 / 64.0) * x60;
                T x82 = x14 * x79;
                T x83 = x18 * x79;
                T x84 = x27 * x79;
                T x85 = x24 * x84;
                T x86 = x16 * (-x25 - x5 * x83 - x85);
                T x87 = x16 * (x35 * x81 + x52 * x81 + x56);
                T x88 = x45 * x83;
                T x89 = x16 * (-x65 * x84 - x75 - x88);
                T x90 = x16 * (x34 * x82 + x45 * x80 + x53);
                T x91 = x63 * x79;
                T x92 = x16 * (x68 * x8 + x72 * x8 + x77);
                T x93 = x16 * (x23 * x91 + x71 * x79 + x73);
                T x94 = x16 * (-x38 * x84 - x49 - x88);
                T x95 = pow(x44, 2);
                T x96 = x58 * x95;
                T x97 = x61 * x95;
                T x98 = x20 * x95;
                T x99 = x24 * x98;
                T x100 = x16 * (-x29 - x46 * x95 - x99);
                T x101 = x16 * (x23 * x97 + x32 * x96 + x69);
                T x102 = x39 * x95;
                T x103 = x16 * (-x102 - x38 * x98 - x42);
                T x104 = x11 * x95;
                T x105 = x16 * (-x102 - x65 * x98 - x66);
                T x106 = x16 * (x104 * x34 + x33 * x95 + x36);
                T x107 = x16 * (-x83 * x95 - x85 - x99);
                H[0] += x16 * (x1 * x6 + x12 * x9 + x15 * x9);
                H[1] += x30;
                H[2] += x37;
                H[3] += x43;
                H[4] += x50;
                H[5] += x54;
                H[6] += x55;
                H[7] += x57;
                H[8] += x30;
                H[9] += x16 * (x3 * x59 + x60 * x62 + x60 * x64);
                H[10] += x67;
                H[11] += x70;
                H[12] += x74;
                H[13] += x76;
                H[14] += x78;
                H[15] += x55;
                H[16] += x37;
                H[17] += x67;
                H[18] += x16 * (x12 * x81 + x5 * x80 + x81 * x82);
                H[19] += x86;
                H[20] += x55;
                H[21] += x87;
                H[22] += x89;
                H[23] += x90;
                H[24] += x43;
                H[25] += x70;
                H[26] += x86;
                H[27] += x16 * (x59 * x79 + x62 * x8 + x8 * x91);
                H[28] += x92;
                H[29] += x55;
                H[30] += x93;
                H[31] += x94;
                H[32] += x50;
                H[33] += x74;
                H[34] += x55;
                H[35] += x92;
                H[36] += x16 * (x3 * x96 + x64 * x8 + x8 * x97);
                H[37] += x100;
                H[38] += x101;
                H[39] += x103;
                H[40] += x54;
                H[41] += x76;
                H[42] += x87;
                H[43] += x55;
                H[44] += x100;
                H[45] += x16 * (x104 * x81 + x15 * x81 + x51 * x95);
                H[46] += x105;
                H[47] += x106;
                H[48] += x55;
                H[49] += x78;
                H[50] += x89;
                H[51] += x93;
                H[52] += x101;
                H[53] += x105;
                H[54] += x16 * (x60 * x91 + x60 * x97 + x79 * x96);
                H[55] += x107;
                H[56] += x57;
                H[57] += x55;
                H[58] += x90;
                H[59] += x94;
                H[60] += x103;
                H[61] += x106;
                H[62] += x107;
                H[63] += x16 * (x104 * x9 + x80 * x95 + x82 * x9);
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
                //	- Result: 8*ADD + 8*ADDAUGMENTEDASSIGNMENT + 32*MUL
                //	- Subexpressions: 12*ADD + 6*DIV + 71*MUL + NEG + 21*SUB
                T x0 = y - 1.0;
                T x1 = z - 1.0;
                T x2 = px[0] - px[6];
                T x3 = (1.0 / 8.0) / x2;
                T x4 = x1 * x3;
                T x5 = u[0] * x0;
                T x6 = y + 1.0;
                T x7 = u[2] * x6;
                T x8 = z + 1.0;
                T x9 = x3 * x8;
                T x10 = u[5] * x0;
                T x11 = u[7] * x6;
                T x12 = 1.0 / (8 * px[0] - 8 * px[6]);
                T x13 = x1 * x12;
                T x14 = u[1] * x0;
                T x15 = u[3] * x6;
                T x16 = x12 * x8;
                T x17 = u[4] * x0;
                T x18 = u[6] * x6;
                T x19 = x10 * x9 + x11 * x9 - x13 * x14 - x13 * x15 - x16 * x17 - x16 * x18 + x4 * x5 + x4 * x7;
                T x20 = x19 * x4;
                T x21 = x - 1.0;
                T x22 = py[0] - py[6];
                T x23 = (1.0 / 8.0) / x22;
                T x24 = x1 * x23;
                T x25 = x21 * x24;
                T x26 = x + 1.0;
                T x27 = x24 * x26;
                T x28 = x23 * x8;
                T x29 = x26 * x28;
                T x30 = x21 * x28;
                T x31 = 1.0 / (8 * py[0] - 8 * py[6]);
                T x32 = x1 * x31;
                T x33 = x26 * x32;
                T x34 = x21 * x32;
                T x35 = x31 * x8;
                T x36 = x21 * x35;
                T x37 = x26 * x35;
                T x38 = u[0] * x25 - u[1] * x33 + u[2] * x27 - u[3] * x34 - u[4] * x36 + u[5] * x29 - u[6] * x37 +
                        u[7] * x30;
                T x39 = pz[0] - pz[6];
                T x40 = (1.0 / 8.0) / x39;
                T x41 = x21 * x40;
                T x42 = x26 * x40;
                T x43 = 1.0 / (8 * pz[0] - 8 * pz[6]);
                T x44 = x26 * x43;
                T x45 = x21 * x43;
                T x46 = x10 * x42 + x11 * x41 - x14 * x44 - x15 * x45 - x17 * x45 - x18 * x44 + x41 * x5 + x42 * x7;
                T x47 = x0 * x46;
                T x48 = -8 * weight * x2 * x22 * x39;
                T x49 = x13 * x19;
                T x50 = x46 * x6;
                T x51 = x16 * x19;
                T x52 = x19 * x9;
                Hx[0] += x48 * (x0 * x20 + x25 * x38 + x41 * x47);
                Hx[1] += x48 * (-x0 * x49 - x33 * x38 - x44 * x47);
                Hx[2] += x48 * (x20 * x6 + x27 * x38 + x42 * x50);
                Hx[3] += x48 * (-x34 * x38 - x45 * x50 - x49 * x6);
                Hx[4] += x48 * (-x0 * x51 - x36 * x38 - x45 * x47);
                Hx[5] += x48 * (x0 * x52 + x29 * x38 + x42 * x47);
                Hx[6] += x48 * (-x37 * x38 - x44 * x50 - x51 * x6);
                Hx[7] += x48 * (x30 * x38 + x41 * x50 + x52 * x6);
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
                //	- Result: 8*ADD + 8*ADDAUGMENTEDASSIGNMENT + 32*MUL
                //	- Subexpressions: 12*ADD + 6*DIV + 71*MUL + NEG + 21*SUB
                T x0 = y - 1.0;
                T x1 = z - 1.0;
                T x2 = px[0] - px[6];
                T x3 = (1.0 / 8.0) / x2;
                T x4 = x1 * x3;
                T x5 = u[0] * x0;
                T x6 = y + 1.0;
                T x7 = u[2] * x6;
                T x8 = z + 1.0;
                T x9 = x3 * x8;
                T x10 = u[5] * x0;
                T x11 = u[7] * x6;
                T x12 = 1.0 / (8 * px[0] - 8 * px[6]);
                T x13 = x1 * x12;
                T x14 = u[1] * x0;
                T x15 = u[3] * x6;
                T x16 = x12 * x8;
                T x17 = u[4] * x0;
                T x18 = u[6] * x6;
                T x19 = x10 * x9 + x11 * x9 - x13 * x14 - x13 * x15 - x16 * x17 - x16 * x18 + x4 * x5 + x4 * x7;
                T x20 = x19 * x4;
                T x21 = x - 1.0;
                T x22 = py[0] - py[6];
                T x23 = (1.0 / 8.0) / x22;
                T x24 = x1 * x23;
                T x25 = x21 * x24;
                T x26 = x + 1.0;
                T x27 = x24 * x26;
                T x28 = x23 * x8;
                T x29 = x26 * x28;
                T x30 = x21 * x28;
                T x31 = 1.0 / (8 * py[0] - 8 * py[6]);
                T x32 = x1 * x31;
                T x33 = x26 * x32;
                T x34 = x21 * x32;
                T x35 = x31 * x8;
                T x36 = x21 * x35;
                T x37 = x26 * x35;
                T x38 = u[0] * x25 - u[1] * x33 + u[2] * x27 - u[3] * x34 - u[4] * x36 + u[5] * x29 - u[6] * x37 +
                        u[7] * x30;
                T x39 = pz[0] - pz[6];
                T x40 = (1.0 / 8.0) / x39;
                T x41 = x21 * x40;
                T x42 = x26 * x40;
                T x43 = 1.0 / (8 * pz[0] - 8 * pz[6]);
                T x44 = x26 * x43;
                T x45 = x21 * x43;
                T x46 = x10 * x42 + x11 * x41 - x14 * x44 - x15 * x45 - x17 * x45 - x18 * x44 + x41 * x5 + x42 * x7;
                T x47 = x0 * x46;
                T x48 = -8 * weight * x2 * x22 * x39;
                T x49 = x13 * x19;
                T x50 = x46 * x6;
                T x51 = x16 * x19;
                T x52 = x19 * x9;
                g[0] += x48 * (x0 * x20 + x25 * x38 + x41 * x47);
                g[1] += x48 * (-x0 * x49 - x33 * x38 - x44 * x47);
                g[2] += x48 * (x20 * x6 + x27 * x38 + x42 * x50);
                g[3] += x48 * (-x34 * x38 - x45 * x50 - x49 * x6);
                g[4] += x48 * (-x0 * x51 - x36 * x38 - x45 * x47);
                g[5] += x48 * (x0 * x52 + x29 * x38 + x42 * x47);
                g[6] += x48 * (-x37 * x38 - x44 * x50 - x51 * x6);
                g[7] += x48 * (x30 * x38 + x41 * x50 + x52 * x6);
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
                //	- Result: 4*ADD + ADDAUGMENTEDASSIGNMENT + 25*MUL + 3*POW
                //	- Subexpressions: 3*ADD + 6*DIV + 29*MUL + 9*SUB
                T x0 = px[0] - px[6];
                T x1 = py[0] - py[6];
                T x2 = pz[0] - pz[6];
                T x3 = z - 1.0;
                T x4 = (1.0 / 8.0) / x0;
                T x5 = x3 * x4;
                T x6 = y - 1.0;
                T x7 = u[0] * x6;
                T x8 = y + 1.0;
                T x9 = u[2] * x8;
                T x10 = z + 1.0;
                T x11 = x10 * x4;
                T x12 = u[5] * x6;
                T x13 = u[7] * x8;
                T x14 = 1.0 / (8 * px[0] - 8 * px[6]);
                T x15 = x14 * x3;
                T x16 = u[1] * x6;
                T x17 = u[3] * x8;
                T x18 = x10 * x14;
                T x19 = u[4] * x6;
                T x20 = u[6] * x8;
                T x21 = x - 1.0;
                T x22 = (1.0 / 8.0) / x1;
                T x23 = x22 * x3;
                T x24 = x + 1.0;
                T x25 = x10 * x22;
                T x26 = 1.0 / (8 * py[0] - 8 * py[6]);
                T x27 = x26 * x3;
                T x28 = x10 * x26;
                T x29 = (1.0 / 8.0) / x2;
                T x30 = x21 * x29;
                T x31 = x24 * x29;
                T x32 = 1.0 / (8 * pz[0] - 8 * pz[6]);
                T x33 = x24 * x32;
                T x34 = x21 * x32;
                e +=
                    -weight * x0 * x1 * x2 *
                    (pow(x11 * x12 + x11 * x13 - x15 * x16 - x15 * x17 - x18 * x19 - x18 * x20 + x5 * x7 + x5 * x9, 2) +
                     pow(x12 * x31 + x13 * x30 - x16 * x33 - x17 * x34 - x19 * x34 - x20 * x33 + x30 * x7 + x31 * x9,
                         2) +
                     pow(u[0] * x21 * x23 - u[1] * x24 * x27 + u[2] * x23 * x24 - u[3] * x21 * x27 - u[4] * x21 * x28 +
                             u[5] * x24 * x25 - u[6] * x24 * x28 + u[7] * x21 * x25,
                         2));
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
                //	- Result: 17*ADD + 73*ADDAUGMENTEDASSIGNMENT + 65*MUL + 3*POW
                //	- Subexpressions: 36*ADD + 15*DIV + 186*MUL + 19*NEG + 12*POW + 47*SUB
                T x0 = px[0] - px[6];
                T x1 = (1.0 / 64.0) / pow(x0, 2);
                T x2 = y - 1.0;
                T x3 = pow(x2, 2);
                T x4 = z - 1.0;
                T x5 = pow(x4, 2);
                T x6 = x3 * x5;
                T x7 = x - 1.0;
                T x8 = pow(x7, 2);
                T x9 = (1.0 / 64.0) * x8;
                T x10 = py[0] - py[6];
                T x11 = pow(x10, -2);
                T x12 = x11 * x5;
                T x13 = pz[0] - pz[6];
                T x14 = pow(x13, -2);
                T x15 = x14 * x3;
                T x16 = -weight * x0 * x10 * x13;
                T x17 = 8 * px[0] - 8 * px[6];
                T x18 = 1.0 / x17;
                T x19 = (1.0 / 8.0) / x0;
                T x20 = x18 * x19;
                T x21 = 8 * py[0] - 8 * py[6];
                T x22 = 1.0 / x21;
                T x23 = 1.0 / x10;
                T x24 = x22 * x23;
                T x25 = x24 * x5;
                T x26 = x + 1.0;
                T x27 = x26 * x7;
                T x28 = (1.0 / 8.0) * x27;
                T x29 = x25 * x28;
                T x30 = 8 * pz[0] - 8 * pz[6];
                T x31 = 1.0 / x30;
                T x32 = 1.0 / x13;
                T x33 = x31 * x32;
                T x34 = x3 * x33;
                T x35 = x28 * x34;
                T x36 = x16 * (-x20 * x6 - x29 - x35);
                T x37 = y + 1.0;
                T x38 = x2 * x37;
                T x39 = x1 * x38;
                T x40 = (1.0 / 64.0) * x27;
                T x41 = x14 * x38;
                T x42 = x40 * x41;
                T x43 = x16 * (x12 * x40 + x39 * x5 + x42);
                T x44 = (1.0 / 8.0) * x8;
                T x45 = x20 * x38;
                T x46 = x45 * x5;
                T x47 = x33 * x38;
                T x48 = x44 * x47;
                T x49 = x16 * (-x25 * x44 - x46 - x48);
                T x50 = z + 1.0;
                T x51 = x4 * x50;
                T x52 = x20 * x3;
                T x53 = x51 * x52;
                T x54 = x24 * x51;
                T x55 = x44 * x54;
                T x56 = x16 * (-x34 * x44 - x53 - x55);
                T x57 = x1 * x3;
                T x58 = x11 * x51;
                T x59 = x40 * x58;
                T x60 = x16 * (x15 * x40 + x51 * x57 + x59);
                T x61 = x16 * (-x28 * x47 - x28 * x54 - x45 * x51);
                T x62 = x39 * x51;
                T x63 = x16 * (x41 * x9 + x58 * x9 + x62);
                T x64 = pow(x17, -2);
                T x65 = x5 * x64;
                T x66 = pow(x26, 2);
                T x67 = pow(x21, -2);
                T x68 = x5 * x67;
                T x69 = pow(x30, -2);
                T x70 = x3 * x69;
                T x71 = (1.0 / 8.0) * x66;
                T x72 = x47 * x71;
                T x73 = x16 * (-x25 * x71 - x46 - x72);
                T x74 = x38 * x69;
                T x75 = x27 * x74;
                T x76 = x16 * (x27 * x68 + x38 * x65 + x75);
                T x77 = x51 * x64;
                T x78 = x51 * x67;
                T x79 = x27 * x78;
                T x80 = x16 * (x27 * x70 + x3 * x77 + x79);
                T x81 = x54 * x71;
                T x82 = x16 * (-x34 * x71 - x53 - x81);
                T x83 = x38 * x77;
                T x84 = x16 * (x66 * x74 + x66 * x78 + x83);
                T x85 = pow(x37, 2);
                T x86 = x1 * x85;
                T x87 = (1.0 / 64.0) * x66;
                T x88 = x14 * x85;
                T x89 = x20 * x85;
                T x90 = x33 * x85;
                T x91 = x28 * x90;
                T x92 = x16 * (-x29 - x5 * x89 - x91);
                T x93 = x16 * (x41 * x87 + x58 * x87 + x62);
                T x94 = x51 * x89;
                T x95 = x16 * (-x71 * x90 - x81 - x94);
                T x96 = x16 * (x40 * x88 + x51 * x86 + x59);
                T x97 = x69 * x85;
                T x98 = x16 * (x74 * x8 + x78 * x8 + x83);
                T x99 = x16 * (x27 * x97 + x77 * x85 + x79);
                T x100 = x16 * (-x44 * x90 - x55 - x94);
                T x101 = pow(x50, 2);
                T x102 = x101 * x64;
                T x103 = x101 * x67;
                T x104 = x101 * x24;
                T x105 = x104 * x28;
                T x106 = x16 * (-x101 * x52 - x105 - x35);
                T x107 = x16 * (x102 * x38 + x103 * x27 + x75);
                T x108 = x101 * x45;
                T x109 = x16 * (-x104 * x44 - x108 - x48);
                T x110 = x101 * x11;
                T x111 = x16 * (-x104 * x71 - x108 - x72);
                T x112 = x16 * (x101 * x39 + x110 * x40 + x42);
                T x113 = x16 * (-x101 * x89 - x105 - x91);
                T x114 = x19 * x4;
                T x115 = u[0] * x2;
                T x116 = u[2] * x37;
                T x117 = x19 * x50;
                T x118 = u[5] * x2;
                T x119 = u[7] * x37;
                T x120 = x18 * x4;
                T x121 = u[1] * x2;
                T x122 = u[3] * x37;
                T x123 = x18 * x50;
                T x124 = u[4] * x2;
                T x125 = u[6] * x37;
                T x126 = x114 * x115 + x114 * x116 + x117 * x118 + x117 * x119 - x120 * x121 - x120 * x122 -
                         x123 * x124 - x123 * x125;
                T x127 = x114 * x126;
                T x128 = (1.0 / 8.0) * x23;
                T x129 = x128 * x4;
                T x130 = x129 * x7;
                T x131 = x129 * x26;
                T x132 = x128 * x50;
                T x133 = x132 * x26;
                T x134 = x132 * x7;
                T x135 = x22 * x4;
                T x136 = x135 * x26;
                T x137 = x135 * x7;
                T x138 = x22 * x50;
                T x139 = x138 * x7;
                T x140 = x138 * x26;
                T x141 = u[0] * x130 - u[1] * x136 + u[2] * x131 - u[3] * x137 - u[4] * x139 + u[5] * x133 -
                         u[6] * x140 + u[7] * x134;
                T x142 = (1.0 / 8.0) * x32;
                T x143 = x142 * x7;
                T x144 = x142 * x26;
                T x145 = x26 * x31;
                T x146 = x31 * x7;
                T x147 = x115 * x143 + x116 * x144 + x118 * x144 + x119 * x143 - x121 * x145 - x122 * x146 -
                         x124 * x146 - x125 * x145;
                T x148 = x147 * x2;
                T x149 = 8 * x16;
                T x150 = x120 * x126;
                T x151 = x147 * x37;
                T x152 = x123 * x126;
                T x153 = x117 * x126;
                H[0] += x16 * (x1 * x6 + x12 * x9 + x15 * x9);
                H[1] += x36;
                H[2] += x43;
                H[3] += x49;
                H[4] += x56;
                H[5] += x60;
                H[6] += x61;
                H[7] += x63;
                H[8] += x36;
                H[9] += x16 * (x3 * x65 + x66 * x68 + x66 * x70);
                H[10] += x73;
                H[11] += x76;
                H[12] += x80;
                H[13] += x82;
                H[14] += x84;
                H[15] += x61;
                H[16] += x43;
                H[17] += x73;
                H[18] += x16 * (x12 * x87 + x5 * x86 + x87 * x88);
                H[19] += x92;
                H[20] += x61;
                H[21] += x93;
                H[22] += x95;
                H[23] += x96;
                H[24] += x49;
                H[25] += x76;
                H[26] += x92;
                H[27] += x16 * (x65 * x85 + x68 * x8 + x8 * x97);
                H[28] += x98;
                H[29] += x61;
                H[30] += x99;
                H[31] += x100;
                H[32] += x56;
                H[33] += x80;
                H[34] += x61;
                H[35] += x98;
                H[36] += x16 * (x102 * x3 + x103 * x8 + x70 * x8);
                H[37] += x106;
                H[38] += x107;
                H[39] += x109;
                H[40] += x60;
                H[41] += x82;
                H[42] += x93;
                H[43] += x61;
                H[44] += x106;
                H[45] += x16 * (x101 * x57 + x110 * x87 + x15 * x87);
                H[46] += x111;
                H[47] += x112;
                H[48] += x61;
                H[49] += x84;
                H[50] += x95;
                H[51] += x99;
                H[52] += x107;
                H[53] += x111;
                H[54] += x16 * (x102 * x85 + x103 * x66 + x66 * x97);
                H[55] += x113;
                H[56] += x63;
                H[57] += x61;
                H[58] += x96;
                H[59] += x100;
                H[60] += x109;
                H[61] += x112;
                H[62] += x113;
                H[63] += x16 * (x101 * x86 + x110 * x9 + x88 * x9);
                g[0] += x149 * (x127 * x2 + x130 * x141 + x143 * x148);
                g[1] += x149 * (-x136 * x141 - x145 * x148 - x150 * x2);
                g[2] += x149 * (x127 * x37 + x131 * x141 + x144 * x151);
                g[3] += x149 * (-x137 * x141 - x146 * x151 - x150 * x37);
                g[4] += x149 * (-x139 * x141 - x146 * x148 - x152 * x2);
                g[5] += x149 * (x133 * x141 + x144 * x148 + x153 * x2);
                g[6] += x149 * (-x140 * x141 - x145 * x151 - x152 * x37);
                g[7] += x149 * (x134 * x141 + x143 * x151 + x153 * x37);
                e += x16 * (pow(x126, 2) + pow(x141, 2) + pow(x147, 2));
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using LaplaceOperatorAxisAlignedHex8 =
            utopia::kokkos::AutoKernel<FE,
                                       utopia::kernels::LaplaceOperator<
                                           utopia::kernels::AxisAlignedHex8<typename FE::Scalar, typename FE::Scalar>>,
                                       3>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_LaplaceOperator_AxisAlignedHex8_3_IMPL_hpp
