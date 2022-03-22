#ifndef UTOPIA_TPL_MATERIAL_LaplaceOperator_Hex8_3_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_LaplaceOperator_Hex8_3_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Hex8_3.hpp"
#include "utopia_material_LaplaceOperator.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of LaplaceOperator for symmetric element pair trial=test=Hex8
         */
        template <typename T, typename GeoT>
        class LaplaceOperator<Hex8<T, GeoT>> {
        public:
            using ElemT = Hex8<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "LaplaceOperator<Hex8>"; }

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
                //	- Result: 8*ADD + 64*ADDAUGMENTEDASSIGNMENT + 8*MUL + 24*POW
                //	- Subexpressions: 118*ADD + 47*DIV + 356*MUL + 5*NEG + 74*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * x3;
                T x5 = x1 * x4;
                T x6 = x0 + 0.5;
                T x7 = x4 * x6;
                T x8 = x2 + 0.5;
                T x9 = (1.0 / 2.0) * x8;
                T x10 = x6 * x9;
                T x11 = x1 * x9;
                T x12 = -px[0] * x5 - px[1] * x7 - px[2] * x10 - px[3] * x11 + (1.0 / 2.0) * px[4] * x1 * x3 +
                        (1.0 / 2.0) * px[5] * x3 * x6 + (1.0 / 2.0) * px[6] * x6 * x8 + (1.0 / 2.0) * px[7] * x1 * x8;
                T x13 = (1.0 / 2.0) * z;
                T x14 = 0.5 - x13;
                T x15 = (1.0 / 2.0) * x14;
                T x16 = x1 * x15;
                T x17 = x15 * x6;
                T x18 = x13 + 0.5;
                T x19 = (1.0 / 2.0) * x18;
                T x20 = x1 * x19;
                T x21 = x19 * x6;
                T x22 = -pz[0] * x16 - pz[1] * x17 + (1.0 / 2.0) * pz[2] * x14 * x6 + (1.0 / 2.0) * pz[3] * x1 * x14 -
                        pz[4] * x20 - pz[5] * x21 + (1.0 / 2.0) * pz[6] * x18 * x6 + (1.0 / 2.0) * pz[7] * x1 * x18;
                T x23 = -px[0] * x16 - px[1] * x17 + (1.0 / 2.0) * px[2] * x14 * x6 + (1.0 / 2.0) * px[3] * x1 * x14 -
                        px[4] * x20 - px[5] * x21 + (1.0 / 2.0) * px[6] * x18 * x6 + (1.0 / 2.0) * px[7] * x1 * x18;
                T x24 = -pz[0] * x5 - pz[1] * x7 - pz[2] * x10 - pz[3] * x11 + (1.0 / 2.0) * pz[4] * x1 * x3 +
                        (1.0 / 2.0) * pz[5] * x3 * x6 + (1.0 / 2.0) * pz[6] * x6 * x8 + (1.0 / 2.0) * pz[7] * x1 * x8;
                T x25 = x12 * x22 - x23 * x24;
                T x26 = (1.0 / 4.0) * x;
                T x27 = x26 - 0.25;
                T x28 = -py[0] * x16 - py[1] * x17 + (1.0 / 2.0) * py[2] * x14 * x6 + (1.0 / 2.0) * py[3] * x1 * x14 -
                        py[4] * x20 - py[5] * x21 + (1.0 / 2.0) * py[6] * x18 * x6 + (1.0 / 2.0) * py[7] * x1 * x18;
                T x29 = x15 * x3;
                T x30 = x15 * x8;
                T x31 = x19 * x3;
                T x32 = x19 * x8;
                T x33 = -pz[0] * x29 + (1.0 / 2.0) * pz[1] * x14 * x3 + (1.0 / 2.0) * pz[2] * x14 * x8 - pz[3] * x30 -
                        pz[4] * x31 + (1.0 / 2.0) * pz[5] * x18 * x3 + (1.0 / 2.0) * pz[6] * x18 * x8 - pz[7] * x32;
                T x34 = x12 * x33;
                T x35 = -py[0] * x29 + (1.0 / 2.0) * py[1] * x14 * x3 + (1.0 / 2.0) * py[2] * x14 * x8 - py[3] * x30 -
                        py[4] * x31 + (1.0 / 2.0) * py[5] * x18 * x3 + (1.0 / 2.0) * py[6] * x18 * x8 - py[7] * x32;
                T x36 = x24 * x35;
                T x37 = -px[0] * x29 + (1.0 / 2.0) * px[1] * x14 * x3 + (1.0 / 2.0) * px[2] * x14 * x8 - px[3] * x30 -
                        px[4] * x31 + (1.0 / 2.0) * px[5] * x18 * x3 + (1.0 / 2.0) * px[6] * x18 * x8 - px[7] * x32;
                T x38 = -py[0] * x5 - py[1] * x7 - py[2] * x10 - py[3] * x11 + (1.0 / 2.0) * py[4] * x1 * x3 +
                        (1.0 / 2.0) * py[5] * x3 * x6 + (1.0 / 2.0) * py[6] * x6 * x8 + (1.0 / 2.0) * py[7] * x1 * x8;
                T x39 = x37 * x38;
                T x40 = x12 * x22 * x35 - x22 * x39 + x23 * x33 * x38 - x23 * x36 + x24 * x28 * x37 - x28 * x34;
                T x41 = 1.0 / x40;
                T x42 = x14 * x41;
                T x43 = x27 * x42;
                T x44 = -x22 * x38 + x24 * x28;
                T x45 = (1.0 / 4.0) * y;
                T x46 = x45 - 0.25;
                T x47 = x42 * x46;
                T x48 = -x12 * x28 + x23 * x38;
                T x49 = x41 * x5;
                T x50 = x25 * x43 + x44 * x47 - x48 * x49;
                T x51 = x24 * x37 - x34;
                T x52 = x33 * x38 - x36;
                T x53 = x12 * x35 - x39;
                T x54 = x43 * x51 + x47 * x52 - x49 * x53;
                T x55 = -x22 * x37 + x23 * x33;
                T x56 = x22 * x35 - x28 * x33;
                T x57 = -x23 * x35 + x28 * x37;
                T x58 = x43 * x55 + x47 * x56 - x49 * x57;
                T x59 = weight * x40;
                T x60 = -x42 * x46;
                T x61 = x26 + 0.25;
                T x62 = -x42 * x61;
                T x63 = x41 * x7;
                T x64 = x25 * x62 + x44 * x60 - x48 * x63;
                T x65 = x51 * x62 + x52 * x60 - x53 * x63;
                T x66 = x55 * x62 + x56 * x60 - x57 * x63;
                T x67 = x59 * (x50 * x64 + x54 * x65 + x58 * x66);
                T x68 = x42 * x61;
                T x69 = x45 + 0.25;
                T x70 = x42 * x69;
                T x71 = x10 * x41;
                T x72 = x48 * x71;
                T x73 = x25 * x68 + x44 * x70 - x72;
                T x74 = x53 * x71;
                T x75 = x51 * x68 + x52 * x70 - x74;
                T x76 = x57 * x71;
                T x77 = x55 * x68 + x56 * x70 - x76;
                T x78 = x59 * (x50 * x73 + x54 * x75 + x58 * x77);
                T x79 = -x27;
                T x80 = x42 * x79;
                T x81 = -x42 * x69;
                T x82 = x11 * x41;
                T x83 = x25 * x80 + x44 * x81 - x48 * x82;
                T x84 = x51 * x80 + x52 * x81 - x53 * x82;
                T x85 = x55 * x80 + x56 * x81 - x57 * x82;
                T x86 = x59 * (x50 * x83 + x54 * x84 + x58 * x85);
                T x87 = x3 * x41;
                T x88 = x79 * x87;
                T x89 = (1.0 / 4.0) * z + 0.25;
                T x90 = -x87 * x89;
                T x91 = x20 * x41;
                T x92 = -x25 * x91 + x44 * x90 + x48 * x88;
                T x93 = -x51 * x91 + x52 * x90 + x53 * x88;
                T x94 = -x55 * x91 + x56 * x90 + x57 * x88;
                T x95 = x59 * (x50 * x92 + x54 * x93 + x58 * x94);
                T x96 = x61 * x87;
                T x97 = x87 * x89;
                T x98 = x21 * x41;
                T x99 = -x25 * x98 + x44 * x97 + x48 * x96;
                T x100 = -x51 * x98 + x52 * x97 + x53 * x96;
                T x101 = -x55 * x98 + x56 * x97 + x57 * x96;
                T x102 = x59 * (x100 * x54 + x101 * x58 + x50 * x99);
                T x103 = x18 * x41;
                T x104 = x103 * x61;
                T x105 = x103 * x69;
                T x106 = x104 * x25 + x105 * x44 + x72;
                T x107 = x104 * x51 + x105 * x52 + x74;
                T x108 = x104 * x55 + x105 * x56 + x76;
                T x109 = x59 * (x106 * x50 + x107 * x54 + x108 * x58);
                T x110 = x1 * x41;
                T x111 = x110 * x69;
                T x112 = x110 * x89;
                T x113 = x32 * x41;
                T x114 = x111 * x48 + x112 * x25 - x113 * x44;
                T x115 = x111 * x53 + x112 * x51 - x113 * x52;
                T x116 = x111 * x57 + x112 * x55 - x113 * x56;
                T x117 = x59 * (x114 * x50 + x115 * x54 + x116 * x58);
                T x118 = x59 * (x64 * x73 + x65 * x75 + x66 * x77);
                T x119 = x59 * (x64 * x83 + x65 * x84 + x66 * x85);
                T x120 = x59 * (x64 * x92 + x65 * x93 + x66 * x94);
                T x121 = x59 * (x100 * x65 + x101 * x66 + x64 * x99);
                T x122 = x59 * (x106 * x64 + x107 * x65 + x108 * x66);
                T x123 = x59 * (x114 * x64 + x115 * x65 + x116 * x66);
                T x124 = x59 * (x73 * x83 + x75 * x84 + x77 * x85);
                T x125 = x59 * (x73 * x92 + x75 * x93 + x77 * x94);
                T x126 = x59 * (x100 * x75 + x101 * x77 + x73 * x99);
                T x127 = x59 * (x106 * x73 + x107 * x75 + x108 * x77);
                T x128 = x59 * (x114 * x73 + x115 * x75 + x116 * x77);
                T x129 = x59 * (x83 * x92 + x84 * x93 + x85 * x94);
                T x130 = x59 * (x100 * x84 + x101 * x85 + x83 * x99);
                T x131 = x59 * (x106 * x83 + x107 * x84 + x108 * x85);
                T x132 = x59 * (x114 * x83 + x115 * x84 + x116 * x85);
                T x133 = x59 * (x100 * x93 + x101 * x94 + x92 * x99);
                T x134 = x59 * (x106 * x92 + x107 * x93 + x108 * x94);
                T x135 = x59 * (x114 * x92 + x115 * x93 + x116 * x94);
                T x136 = x59 * (x100 * x107 + x101 * x108 + x106 * x99);
                T x137 = x59 * (x100 * x115 + x101 * x116 + x114 * x99);
                T x138 = x59 * (x106 * x114 + x107 * x115 + x108 * x116);
                H[0] += x59 * (pow(x50, 2) + pow(x54, 2) + pow(x58, 2));
                H[1] += x67;
                H[2] += x78;
                H[3] += x86;
                H[4] += x95;
                H[5] += x102;
                H[6] += x109;
                H[7] += x117;
                H[8] += x67;
                H[9] += x59 * (pow(x64, 2) + pow(x65, 2) + pow(x66, 2));
                H[10] += x118;
                H[11] += x119;
                H[12] += x120;
                H[13] += x121;
                H[14] += x122;
                H[15] += x123;
                H[16] += x78;
                H[17] += x118;
                H[18] += x59 * (pow(x73, 2) + pow(x75, 2) + pow(x77, 2));
                H[19] += x124;
                H[20] += x125;
                H[21] += x126;
                H[22] += x127;
                H[23] += x128;
                H[24] += x86;
                H[25] += x119;
                H[26] += x124;
                H[27] += x59 * (pow(x83, 2) + pow(x84, 2) + pow(x85, 2));
                H[28] += x129;
                H[29] += x130;
                H[30] += x131;
                H[31] += x132;
                H[32] += x95;
                H[33] += x120;
                H[34] += x125;
                H[35] += x129;
                H[36] += x59 * (pow(x92, 2) + pow(x93, 2) + pow(x94, 2));
                H[37] += x133;
                H[38] += x134;
                H[39] += x135;
                H[40] += x102;
                H[41] += x121;
                H[42] += x126;
                H[43] += x130;
                H[44] += x133;
                H[45] += x59 * (pow(x100, 2) + pow(x101, 2) + pow(x99, 2));
                H[46] += x136;
                H[47] += x137;
                H[48] += x109;
                H[49] += x122;
                H[50] += x127;
                H[51] += x131;
                H[52] += x134;
                H[53] += x136;
                H[54] += x59 * (pow(x106, 2) + pow(x107, 2) + pow(x108, 2));
                H[55] += x138;
                H[56] += x117;
                H[57] += x123;
                H[58] += x128;
                H[59] += x132;
                H[60] += x135;
                H[61] += x137;
                H[62] += x138;
                H[63] += x59 * (pow(x114, 2) + pow(x115, 2) + pow(x116, 2));
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
                //	- Subexpressions: 83*ADD + 47*DIV + 269*MUL + 5*NEG + 74*SUB
                T x0 = (1.0 / 4.0) * x;
                T x1 = x0 - 0.25;
                T x2 = (1.0 / 2.0) * x;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * y;
                T x5 = 0.5 - x4;
                T x6 = (1.0 / 2.0) * x5;
                T x7 = x3 * x6;
                T x8 = x2 + 0.5;
                T x9 = x6 * x8;
                T x10 = x4 + 0.5;
                T x11 = (1.0 / 2.0) * x10;
                T x12 = x11 * x8;
                T x13 = x11 * x3;
                T x14 = -px[0] * x7 - px[1] * x9 - px[2] * x12 - px[3] * x13 + (1.0 / 2.0) * px[4] * x3 * x5 +
                        (1.0 / 2.0) * px[5] * x5 * x8 + (1.0 / 2.0) * px[6] * x10 * x8 + (1.0 / 2.0) * px[7] * x10 * x3;
                T x15 = (1.0 / 2.0) * z;
                T x16 = 0.5 - x15;
                T x17 = (1.0 / 2.0) * x16;
                T x18 = x17 * x3;
                T x19 = x17 * x8;
                T x20 = x15 + 0.5;
                T x21 = (1.0 / 2.0) * x20;
                T x22 = x21 * x3;
                T x23 = x21 * x8;
                T x24 = -pz[0] * x18 - pz[1] * x19 + (1.0 / 2.0) * pz[2] * x16 * x8 + (1.0 / 2.0) * pz[3] * x16 * x3 -
                        pz[4] * x22 - pz[5] * x23 + (1.0 / 2.0) * pz[6] * x20 * x8 + (1.0 / 2.0) * pz[7] * x20 * x3;
                T x25 = -px[0] * x18 - px[1] * x19 + (1.0 / 2.0) * px[2] * x16 * x8 + (1.0 / 2.0) * px[3] * x16 * x3 -
                        px[4] * x22 - px[5] * x23 + (1.0 / 2.0) * px[6] * x20 * x8 + (1.0 / 2.0) * px[7] * x20 * x3;
                T x26 = -pz[0] * x7 - pz[1] * x9 - pz[2] * x12 - pz[3] * x13 + (1.0 / 2.0) * pz[4] * x3 * x5 +
                        (1.0 / 2.0) * pz[5] * x5 * x8 + (1.0 / 2.0) * pz[6] * x10 * x8 + (1.0 / 2.0) * pz[7] * x10 * x3;
                T x27 = x14 * x24 - x25 * x26;
                T x28 = -py[0] * x18 - py[1] * x19 + (1.0 / 2.0) * py[2] * x16 * x8 + (1.0 / 2.0) * py[3] * x16 * x3 -
                        py[4] * x22 - py[5] * x23 + (1.0 / 2.0) * py[6] * x20 * x8 + (1.0 / 2.0) * py[7] * x20 * x3;
                T x29 = x17 * x5;
                T x30 = x10 * x17;
                T x31 = x21 * x5;
                T x32 = x10 * x21;
                T x33 = -pz[0] * x29 + (1.0 / 2.0) * pz[1] * x16 * x5 + (1.0 / 2.0) * pz[2] * x10 * x16 - pz[3] * x30 -
                        pz[4] * x31 + (1.0 / 2.0) * pz[5] * x20 * x5 + (1.0 / 2.0) * pz[6] * x10 * x20 - pz[7] * x32;
                T x34 = x14 * x33;
                T x35 = -py[0] * x29 + (1.0 / 2.0) * py[1] * x16 * x5 + (1.0 / 2.0) * py[2] * x10 * x16 - py[3] * x30 -
                        py[4] * x31 + (1.0 / 2.0) * py[5] * x20 * x5 + (1.0 / 2.0) * py[6] * x10 * x20 - py[7] * x32;
                T x36 = x26 * x35;
                T x37 = -px[0] * x29 + (1.0 / 2.0) * px[1] * x16 * x5 + (1.0 / 2.0) * px[2] * x10 * x16 - px[3] * x30 -
                        px[4] * x31 + (1.0 / 2.0) * px[5] * x20 * x5 + (1.0 / 2.0) * px[6] * x10 * x20 - px[7] * x32;
                T x38 = -py[0] * x7 - py[1] * x9 - py[2] * x12 - py[3] * x13 + (1.0 / 2.0) * py[4] * x3 * x5 +
                        (1.0 / 2.0) * py[5] * x5 * x8 + (1.0 / 2.0) * py[6] * x10 * x8 + (1.0 / 2.0) * py[7] * x10 * x3;
                T x39 = x37 * x38;
                T x40 = x14 * x24 * x35 - x24 * x39 + x25 * x33 * x38 - x25 * x36 + x26 * x28 * x37 - x28 * x34;
                T x41 = 1.0 / x40;
                T x42 = x16 * x41;
                T x43 = x27 * x42;
                T x44 = (1.0 / 4.0) * y;
                T x45 = x44 - 0.25;
                T x46 = -x24 * x38 + x26 * x28;
                T x47 = x42 * x46;
                T x48 = -x14 * x28 + x25 * x38;
                T x49 = x41 * x48;
                T x50 = x1 * x43 + x45 * x47 - x49 * x7;
                T x51 = -x45;
                T x52 = x0 + 0.25;
                T x53 = -x52;
                T x54 = x43 * x53 + x47 * x51 - x49 * x9;
                T x55 = x44 + 0.25;
                T x56 = x12 * x49;
                T x57 = x43 * x52 + x47 * x55 - x56;
                T x58 = -x1;
                T x59 = -x55;
                T x60 = -x13 * x49 + x43 * x58 + x47 * x59;
                T x61 = x41 * x5;
                T x62 = x48 * x61;
                T x63 = (1.0 / 4.0) * z + 0.25;
                T x64 = -x63;
                T x65 = x46 * x61;
                T x66 = x27 * x41;
                T x67 = -x22 * x66 + x58 * x62 + x64 * x65;
                T x68 = -x23 * x66 + x52 * x62 + x63 * x65;
                T x69 = x20 * x41;
                T x70 = x52 * x69;
                T x71 = x55 * x69;
                T x72 = x27 * x70 + x46 * x71 + x56;
                T x73 = x3 * x41;
                T x74 = x55 * x73;
                T x75 = x63 * x73;
                T x76 = x32 * x41;
                T x77 = x27 * x75 - x46 * x76 + x48 * x74;
                T x78 = u[0] * x50 + u[1] * x54 + u[2] * x57 + u[3] * x60 + u[4] * x67 + u[5] * x68 + u[6] * x72 +
                        u[7] * x77;
                T x79 = x26 * x37 - x34;
                T x80 = x42 * x79;
                T x81 = x33 * x38 - x36;
                T x82 = x42 * x81;
                T x83 = x14 * x35 - x39;
                T x84 = x41 * x83;
                T x85 = x1 * x80 + x45 * x82 - x7 * x84;
                T x86 = x51 * x82 + x53 * x80 - x84 * x9;
                T x87 = x12 * x84;
                T x88 = x52 * x80 + x55 * x82 - x87;
                T x89 = -x13 * x84 + x58 * x80 + x59 * x82;
                T x90 = x61 * x83;
                T x91 = x61 * x81;
                T x92 = x41 * x79;
                T x93 = -x22 * x92 + x58 * x90 + x64 * x91;
                T x94 = -x23 * x92 + x52 * x90 + x63 * x91;
                T x95 = x70 * x79 + x71 * x81 + x87;
                T x96 = x74 * x83 + x75 * x79 - x76 * x81;
                T x97 = u[0] * x85 + u[1] * x86 + u[2] * x88 + u[3] * x89 + u[4] * x93 + u[5] * x94 + u[6] * x95 +
                        u[7] * x96;
                T x98 = -x24 * x37 + x25 * x33;
                T x99 = x42 * x98;
                T x100 = x24 * x35 - x28 * x33;
                T x101 = x100 * x42;
                T x102 = -x25 * x35 + x28 * x37;
                T x103 = x102 * x41;
                T x104 = x1 * x99 + x101 * x45 - x103 * x7;
                T x105 = x101 * x51 - x103 * x9 + x53 * x99;
                T x106 = x103 * x12;
                T x107 = x101 * x55 - x106 + x52 * x99;
                T x108 = x101 * x59 - x103 * x13 + x58 * x99;
                T x109 = x102 * x61;
                T x110 = x100 * x61;
                T x111 = x41 * x98;
                T x112 = x109 * x58 + x110 * x64 - x111 * x22;
                T x113 = x109 * x52 + x110 * x63 - x111 * x23;
                T x114 = x100 * x71 + x106 + x70 * x98;
                T x115 = -x100 * x76 + x102 * x74 + x75 * x98;
                T x116 = u[0] * x104 + u[1] * x105 + u[2] * x107 + u[3] * x108 + u[4] * x112 + u[5] * x113 +
                         u[6] * x114 + u[7] * x115;
                T x117 = 8 * weight * x40;
                Hx[0] += x117 * (x104 * x116 + x50 * x78 + x85 * x97);
                Hx[1] += x117 * (x105 * x116 + x54 * x78 + x86 * x97);
                Hx[2] += x117 * (x107 * x116 + x57 * x78 + x88 * x97);
                Hx[3] += x117 * (x108 * x116 + x60 * x78 + x89 * x97);
                Hx[4] += x117 * (x112 * x116 + x67 * x78 + x93 * x97);
                Hx[5] += x117 * (x113 * x116 + x68 * x78 + x94 * x97);
                Hx[6] += x117 * (x114 * x116 + x72 * x78 + x95 * x97);
                Hx[7] += x117 * (x115 * x116 + x77 * x78 + x96 * x97);
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
                //	- Subexpressions: 83*ADD + 47*DIV + 269*MUL + 5*NEG + 74*SUB
                T x0 = (1.0 / 4.0) * x;
                T x1 = x0 - 0.25;
                T x2 = (1.0 / 2.0) * x;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * y;
                T x5 = 0.5 - x4;
                T x6 = (1.0 / 2.0) * x5;
                T x7 = x3 * x6;
                T x8 = x2 + 0.5;
                T x9 = x6 * x8;
                T x10 = x4 + 0.5;
                T x11 = (1.0 / 2.0) * x10;
                T x12 = x11 * x8;
                T x13 = x11 * x3;
                T x14 = -px[0] * x7 - px[1] * x9 - px[2] * x12 - px[3] * x13 + (1.0 / 2.0) * px[4] * x3 * x5 +
                        (1.0 / 2.0) * px[5] * x5 * x8 + (1.0 / 2.0) * px[6] * x10 * x8 + (1.0 / 2.0) * px[7] * x10 * x3;
                T x15 = (1.0 / 2.0) * z;
                T x16 = 0.5 - x15;
                T x17 = (1.0 / 2.0) * x16;
                T x18 = x17 * x3;
                T x19 = x17 * x8;
                T x20 = x15 + 0.5;
                T x21 = (1.0 / 2.0) * x20;
                T x22 = x21 * x3;
                T x23 = x21 * x8;
                T x24 = -pz[0] * x18 - pz[1] * x19 + (1.0 / 2.0) * pz[2] * x16 * x8 + (1.0 / 2.0) * pz[3] * x16 * x3 -
                        pz[4] * x22 - pz[5] * x23 + (1.0 / 2.0) * pz[6] * x20 * x8 + (1.0 / 2.0) * pz[7] * x20 * x3;
                T x25 = -px[0] * x18 - px[1] * x19 + (1.0 / 2.0) * px[2] * x16 * x8 + (1.0 / 2.0) * px[3] * x16 * x3 -
                        px[4] * x22 - px[5] * x23 + (1.0 / 2.0) * px[6] * x20 * x8 + (1.0 / 2.0) * px[7] * x20 * x3;
                T x26 = -pz[0] * x7 - pz[1] * x9 - pz[2] * x12 - pz[3] * x13 + (1.0 / 2.0) * pz[4] * x3 * x5 +
                        (1.0 / 2.0) * pz[5] * x5 * x8 + (1.0 / 2.0) * pz[6] * x10 * x8 + (1.0 / 2.0) * pz[7] * x10 * x3;
                T x27 = x14 * x24 - x25 * x26;
                T x28 = -py[0] * x18 - py[1] * x19 + (1.0 / 2.0) * py[2] * x16 * x8 + (1.0 / 2.0) * py[3] * x16 * x3 -
                        py[4] * x22 - py[5] * x23 + (1.0 / 2.0) * py[6] * x20 * x8 + (1.0 / 2.0) * py[7] * x20 * x3;
                T x29 = x17 * x5;
                T x30 = x10 * x17;
                T x31 = x21 * x5;
                T x32 = x10 * x21;
                T x33 = -pz[0] * x29 + (1.0 / 2.0) * pz[1] * x16 * x5 + (1.0 / 2.0) * pz[2] * x10 * x16 - pz[3] * x30 -
                        pz[4] * x31 + (1.0 / 2.0) * pz[5] * x20 * x5 + (1.0 / 2.0) * pz[6] * x10 * x20 - pz[7] * x32;
                T x34 = x14 * x33;
                T x35 = -py[0] * x29 + (1.0 / 2.0) * py[1] * x16 * x5 + (1.0 / 2.0) * py[2] * x10 * x16 - py[3] * x30 -
                        py[4] * x31 + (1.0 / 2.0) * py[5] * x20 * x5 + (1.0 / 2.0) * py[6] * x10 * x20 - py[7] * x32;
                T x36 = x26 * x35;
                T x37 = -px[0] * x29 + (1.0 / 2.0) * px[1] * x16 * x5 + (1.0 / 2.0) * px[2] * x10 * x16 - px[3] * x30 -
                        px[4] * x31 + (1.0 / 2.0) * px[5] * x20 * x5 + (1.0 / 2.0) * px[6] * x10 * x20 - px[7] * x32;
                T x38 = -py[0] * x7 - py[1] * x9 - py[2] * x12 - py[3] * x13 + (1.0 / 2.0) * py[4] * x3 * x5 +
                        (1.0 / 2.0) * py[5] * x5 * x8 + (1.0 / 2.0) * py[6] * x10 * x8 + (1.0 / 2.0) * py[7] * x10 * x3;
                T x39 = x37 * x38;
                T x40 = x14 * x24 * x35 - x24 * x39 + x25 * x33 * x38 - x25 * x36 + x26 * x28 * x37 - x28 * x34;
                T x41 = 1.0 / x40;
                T x42 = x16 * x41;
                T x43 = x27 * x42;
                T x44 = (1.0 / 4.0) * y;
                T x45 = x44 - 0.25;
                T x46 = -x24 * x38 + x26 * x28;
                T x47 = x42 * x46;
                T x48 = -x14 * x28 + x25 * x38;
                T x49 = x41 * x48;
                T x50 = x1 * x43 + x45 * x47 - x49 * x7;
                T x51 = -x45;
                T x52 = x0 + 0.25;
                T x53 = -x52;
                T x54 = x43 * x53 + x47 * x51 - x49 * x9;
                T x55 = x44 + 0.25;
                T x56 = x12 * x49;
                T x57 = x43 * x52 + x47 * x55 - x56;
                T x58 = -x1;
                T x59 = -x55;
                T x60 = -x13 * x49 + x43 * x58 + x47 * x59;
                T x61 = x41 * x5;
                T x62 = x48 * x61;
                T x63 = (1.0 / 4.0) * z + 0.25;
                T x64 = -x63;
                T x65 = x46 * x61;
                T x66 = x27 * x41;
                T x67 = -x22 * x66 + x58 * x62 + x64 * x65;
                T x68 = -x23 * x66 + x52 * x62 + x63 * x65;
                T x69 = x20 * x41;
                T x70 = x52 * x69;
                T x71 = x55 * x69;
                T x72 = x27 * x70 + x46 * x71 + x56;
                T x73 = x3 * x41;
                T x74 = x55 * x73;
                T x75 = x63 * x73;
                T x76 = x32 * x41;
                T x77 = x27 * x75 - x46 * x76 + x48 * x74;
                T x78 = u[0] * x50 + u[1] * x54 + u[2] * x57 + u[3] * x60 + u[4] * x67 + u[5] * x68 + u[6] * x72 +
                        u[7] * x77;
                T x79 = x26 * x37 - x34;
                T x80 = x42 * x79;
                T x81 = x33 * x38 - x36;
                T x82 = x42 * x81;
                T x83 = x14 * x35 - x39;
                T x84 = x41 * x83;
                T x85 = x1 * x80 + x45 * x82 - x7 * x84;
                T x86 = x51 * x82 + x53 * x80 - x84 * x9;
                T x87 = x12 * x84;
                T x88 = x52 * x80 + x55 * x82 - x87;
                T x89 = -x13 * x84 + x58 * x80 + x59 * x82;
                T x90 = x61 * x83;
                T x91 = x61 * x81;
                T x92 = x41 * x79;
                T x93 = -x22 * x92 + x58 * x90 + x64 * x91;
                T x94 = -x23 * x92 + x52 * x90 + x63 * x91;
                T x95 = x70 * x79 + x71 * x81 + x87;
                T x96 = x74 * x83 + x75 * x79 - x76 * x81;
                T x97 = u[0] * x85 + u[1] * x86 + u[2] * x88 + u[3] * x89 + u[4] * x93 + u[5] * x94 + u[6] * x95 +
                        u[7] * x96;
                T x98 = -x24 * x37 + x25 * x33;
                T x99 = x42 * x98;
                T x100 = x24 * x35 - x28 * x33;
                T x101 = x100 * x42;
                T x102 = -x25 * x35 + x28 * x37;
                T x103 = x102 * x41;
                T x104 = x1 * x99 + x101 * x45 - x103 * x7;
                T x105 = x101 * x51 - x103 * x9 + x53 * x99;
                T x106 = x103 * x12;
                T x107 = x101 * x55 - x106 + x52 * x99;
                T x108 = x101 * x59 - x103 * x13 + x58 * x99;
                T x109 = x102 * x61;
                T x110 = x100 * x61;
                T x111 = x41 * x98;
                T x112 = x109 * x58 + x110 * x64 - x111 * x22;
                T x113 = x109 * x52 + x110 * x63 - x111 * x23;
                T x114 = x100 * x71 + x106 + x70 * x98;
                T x115 = -x100 * x76 + x102 * x74 + x75 * x98;
                T x116 = u[0] * x104 + u[1] * x105 + u[2] * x107 + u[3] * x108 + u[4] * x112 + u[5] * x113 +
                         u[6] * x114 + u[7] * x115;
                T x117 = 8 * weight * x40;
                g[0] += x117 * (x104 * x116 + x50 * x78 + x85 * x97);
                g[1] += x117 * (x105 * x116 + x54 * x78 + x86 * x97);
                g[2] += x117 * (x107 * x116 + x57 * x78 + x88 * x97);
                g[3] += x117 * (x108 * x116 + x60 * x78 + x89 * x97);
                g[4] += x117 * (x112 * x116 + x67 * x78 + x93 * x97);
                g[5] += x117 * (x113 * x116 + x68 * x78 + x94 * x97);
                g[6] += x117 * (x114 * x116 + x72 * x78 + x95 * x97);
                g[7] += x117 * (x115 * x116 + x77 * x78 + x96 * x97);
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
                //	- Result: 28*ADD + ADDAUGMENTEDASSIGNMENT + 94*MUL + 3*POW
                //	- Subexpressions: 35*ADD + 47*DIV + 177*MUL + 5*NEG + 53*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * z;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * x3;
                T x5 = x1 * x4;
                T x6 = x0 + 0.5;
                T x7 = x4 * x6;
                T x8 = x2 + 0.5;
                T x9 = (1.0 / 2.0) * x8;
                T x10 = x1 * x9;
                T x11 = x6 * x9;
                T x12 = -py[0] * x5 - py[1] * x7 + (1.0 / 2.0) * py[2] * x3 * x6 + (1.0 / 2.0) * py[3] * x1 * x3 -
                        py[4] * x10 - py[5] * x11 + (1.0 / 2.0) * py[6] * x6 * x8 + (1.0 / 2.0) * py[7] * x1 * x8;
                T x13 = (1.0 / 2.0) * y;
                T x14 = 0.5 - x13;
                T x15 = (1.0 / 2.0) * x14;
                T x16 = x1 * x15;
                T x17 = x15 * x6;
                T x18 = x13 + 0.5;
                T x19 = (1.0 / 2.0) * x18;
                T x20 = x19 * x6;
                T x21 = x1 * x19;
                T x22 = -px[0] * x16 - px[1] * x17 - px[2] * x20 - px[3] * x21 + (1.0 / 2.0) * px[4] * x1 * x14 +
                        (1.0 / 2.0) * px[5] * x14 * x6 + (1.0 / 2.0) * px[6] * x18 * x6 +
                        (1.0 / 2.0) * px[7] * x1 * x18;
                T x23 = x14 * x4;
                T x24 = x18 * x4;
                T x25 = x14 * x9;
                T x26 = x18 * x9;
                T x27 = -pz[0] * x23 + (1.0 / 2.0) * pz[1] * x14 * x3 + (1.0 / 2.0) * pz[2] * x18 * x3 - pz[3] * x24 -
                        pz[4] * x25 + (1.0 / 2.0) * pz[5] * x14 * x8 + (1.0 / 2.0) * pz[6] * x18 * x8 - pz[7] * x26;
                T x28 = x22 * x27;
                T x29 = -px[0] * x5 - px[1] * x7 + (1.0 / 2.0) * px[2] * x3 * x6 + (1.0 / 2.0) * px[3] * x1 * x3 -
                        px[4] * x10 - px[5] * x11 + (1.0 / 2.0) * px[6] * x6 * x8 + (1.0 / 2.0) * px[7] * x1 * x8;
                T x30 = -py[0] * x23 + (1.0 / 2.0) * py[1] * x14 * x3 + (1.0 / 2.0) * py[2] * x18 * x3 - py[3] * x24 -
                        py[4] * x25 + (1.0 / 2.0) * py[5] * x14 * x8 + (1.0 / 2.0) * py[6] * x18 * x8 - py[7] * x26;
                T x31 = -pz[0] * x16 - pz[1] * x17 - pz[2] * x20 - pz[3] * x21 + (1.0 / 2.0) * pz[4] * x1 * x14 +
                        (1.0 / 2.0) * pz[5] * x14 * x6 + (1.0 / 2.0) * pz[6] * x18 * x6 +
                        (1.0 / 2.0) * pz[7] * x1 * x18;
                T x32 = x30 * x31;
                T x33 = -pz[0] * x5 - pz[1] * x7 + (1.0 / 2.0) * pz[2] * x3 * x6 + (1.0 / 2.0) * pz[3] * x1 * x3 -
                        pz[4] * x10 - pz[5] * x11 + (1.0 / 2.0) * pz[6] * x6 * x8 + (1.0 / 2.0) * pz[7] * x1 * x8;
                T x34 = -px[0] * x23 + (1.0 / 2.0) * px[1] * x14 * x3 + (1.0 / 2.0) * px[2] * x18 * x3 - px[3] * x24 -
                        px[4] * x25 + (1.0 / 2.0) * px[5] * x14 * x8 + (1.0 / 2.0) * px[6] * x18 * x8 - px[7] * x26;
                T x35 = -py[0] * x16 - py[1] * x17 - py[2] * x20 - py[3] * x21 + (1.0 / 2.0) * py[4] * x1 * x14 +
                        (1.0 / 2.0) * py[5] * x14 * x6 + (1.0 / 2.0) * py[6] * x18 * x6 +
                        (1.0 / 2.0) * py[7] * x1 * x18;
                T x36 = x34 * x35;
                T x37 = -x12 * x28 + x12 * x31 * x34 + x22 * x30 * x33 + x27 * x29 * x35 - x29 * x32 - x33 * x36;
                T x38 = (1.0 / 4.0) * x;
                T x39 = x38 - 0.25;
                T x40 = x22 * x33 - x29 * x31;
                T x41 = 1.0 / x37;
                T x42 = x3 * x41;
                T x43 = x40 * x42;
                T x44 = (1.0 / 4.0) * y;
                T x45 = x44 - 0.25;
                T x46 = x12 * x31 - x33 * x35;
                T x47 = x42 * x46;
                T x48 = -x12 * x22 + x29 * x35;
                T x49 = x41 * x48;
                T x50 = -x45;
                T x51 = x38 + 0.25;
                T x52 = -x51;
                T x53 = x44 + 0.25;
                T x54 = x20 * x49;
                T x55 = -x39;
                T x56 = -x53;
                T x57 = x14 * x41;
                T x58 = x48 * x57;
                T x59 = (1.0 / 4.0) * z + 0.25;
                T x60 = -x59;
                T x61 = x46 * x57;
                T x62 = x40 * x41;
                T x63 = x41 * x8;
                T x64 = x51 * x63;
                T x65 = x53 * x63;
                T x66 = x1 * x41;
                T x67 = x53 * x66;
                T x68 = x59 * x66;
                T x69 = x26 * x41;
                T x70 = -x28 + x31 * x34;
                T x71 = x42 * x70;
                T x72 = x27 * x35 - x32;
                T x73 = x42 * x72;
                T x74 = x22 * x30 - x36;
                T x75 = x41 * x74;
                T x76 = x20 * x75;
                T x77 = x57 * x74;
                T x78 = x57 * x72;
                T x79 = x41 * x70;
                T x80 = x27 * x29 - x33 * x34;
                T x81 = x42 * x80;
                T x82 = -x12 * x27 + x30 * x33;
                T x83 = x42 * x82;
                T x84 = x12 * x34 - x29 * x30;
                T x85 = x41 * x84;
                T x86 = x20 * x85;
                T x87 = x57 * x84;
                T x88 = x57 * x82;
                T x89 = x41 * x80;
                e +=
                    weight * x37 *
                    (pow(u[0] * (-x16 * x49 + x39 * x43 + x45 * x47) + u[1] * (-x17 * x49 + x43 * x52 + x47 * x50) +
                             u[2] * (x43 * x51 + x47 * x53 - x54) + u[3] * (-x21 * x49 + x43 * x55 + x47 * x56) +
                             u[4] * (-x10 * x62 + x55 * x58 + x60 * x61) + u[5] * (-x11 * x62 + x51 * x58 + x59 * x61) +
                             u[6] * (x40 * x64 + x46 * x65 + x54) + u[7] * (x40 * x68 - x46 * x69 + x48 * x67),
                         2) +
                     pow(u[0] * (-x16 * x75 + x39 * x71 + x45 * x73) + u[1] * (-x17 * x75 + x50 * x73 + x52 * x71) +
                             u[2] * (x51 * x71 + x53 * x73 - x76) + u[3] * (-x21 * x75 + x55 * x71 + x56 * x73) +
                             u[4] * (-x10 * x79 + x55 * x77 + x60 * x78) + u[5] * (-x11 * x79 + x51 * x77 + x59 * x78) +
                             u[6] * (x64 * x70 + x65 * x72 + x76) + u[7] * (x67 * x74 + x68 * x70 - x69 * x72),
                         2) +
                     pow(u[0] * (-x16 * x85 + x39 * x81 + x45 * x83) + u[1] * (-x17 * x85 + x50 * x83 + x52 * x81) +
                             u[2] * (x51 * x81 + x53 * x83 - x86) + u[3] * (-x21 * x85 + x55 * x81 + x56 * x83) +
                             u[4] * (-x10 * x89 + x55 * x87 + x60 * x88) + u[5] * (-x11 * x89 + x51 * x87 + x59 * x88) +
                             u[6] * (x64 * x80 + x65 * x82 + x86) + u[7] * (x67 * x84 + x68 * x80 - x69 * x82),
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
                //	- Result: 17*ADD + 73*ADDAUGMENTEDASSIGNMENT + 41*MUL + 27*POW
                //	- Subexpressions: 139*ADD + 47*DIV + 381*MUL + 5*NEG + 74*SUB
                T x0 = (1.0 / 2.0) * x;
                T x1 = 0.5 - x0;
                T x2 = (1.0 / 2.0) * y;
                T x3 = 0.5 - x2;
                T x4 = (1.0 / 2.0) * x3;
                T x5 = x1 * x4;
                T x6 = x0 + 0.5;
                T x7 = x4 * x6;
                T x8 = x2 + 0.5;
                T x9 = (1.0 / 2.0) * x8;
                T x10 = x6 * x9;
                T x11 = x1 * x9;
                T x12 = -px[0] * x5 - px[1] * x7 - px[2] * x10 - px[3] * x11 + (1.0 / 2.0) * px[4] * x1 * x3 +
                        (1.0 / 2.0) * px[5] * x3 * x6 + (1.0 / 2.0) * px[6] * x6 * x8 + (1.0 / 2.0) * px[7] * x1 * x8;
                T x13 = (1.0 / 2.0) * z;
                T x14 = 0.5 - x13;
                T x15 = (1.0 / 2.0) * x14;
                T x16 = x1 * x15;
                T x17 = x15 * x6;
                T x18 = x13 + 0.5;
                T x19 = (1.0 / 2.0) * x18;
                T x20 = x1 * x19;
                T x21 = x19 * x6;
                T x22 = -pz[0] * x16 - pz[1] * x17 + (1.0 / 2.0) * pz[2] * x14 * x6 + (1.0 / 2.0) * pz[3] * x1 * x14 -
                        pz[4] * x20 - pz[5] * x21 + (1.0 / 2.0) * pz[6] * x18 * x6 + (1.0 / 2.0) * pz[7] * x1 * x18;
                T x23 = -px[0] * x16 - px[1] * x17 + (1.0 / 2.0) * px[2] * x14 * x6 + (1.0 / 2.0) * px[3] * x1 * x14 -
                        px[4] * x20 - px[5] * x21 + (1.0 / 2.0) * px[6] * x18 * x6 + (1.0 / 2.0) * px[7] * x1 * x18;
                T x24 = -pz[0] * x5 - pz[1] * x7 - pz[2] * x10 - pz[3] * x11 + (1.0 / 2.0) * pz[4] * x1 * x3 +
                        (1.0 / 2.0) * pz[5] * x3 * x6 + (1.0 / 2.0) * pz[6] * x6 * x8 + (1.0 / 2.0) * pz[7] * x1 * x8;
                T x25 = x12 * x22 - x23 * x24;
                T x26 = (1.0 / 4.0) * x;
                T x27 = x26 - 0.25;
                T x28 = -py[0] * x16 - py[1] * x17 + (1.0 / 2.0) * py[2] * x14 * x6 + (1.0 / 2.0) * py[3] * x1 * x14 -
                        py[4] * x20 - py[5] * x21 + (1.0 / 2.0) * py[6] * x18 * x6 + (1.0 / 2.0) * py[7] * x1 * x18;
                T x29 = x15 * x3;
                T x30 = x15 * x8;
                T x31 = x19 * x3;
                T x32 = x19 * x8;
                T x33 = -pz[0] * x29 + (1.0 / 2.0) * pz[1] * x14 * x3 + (1.0 / 2.0) * pz[2] * x14 * x8 - pz[3] * x30 -
                        pz[4] * x31 + (1.0 / 2.0) * pz[5] * x18 * x3 + (1.0 / 2.0) * pz[6] * x18 * x8 - pz[7] * x32;
                T x34 = x12 * x33;
                T x35 = -py[0] * x29 + (1.0 / 2.0) * py[1] * x14 * x3 + (1.0 / 2.0) * py[2] * x14 * x8 - py[3] * x30 -
                        py[4] * x31 + (1.0 / 2.0) * py[5] * x18 * x3 + (1.0 / 2.0) * py[6] * x18 * x8 - py[7] * x32;
                T x36 = x24 * x35;
                T x37 = -px[0] * x29 + (1.0 / 2.0) * px[1] * x14 * x3 + (1.0 / 2.0) * px[2] * x14 * x8 - px[3] * x30 -
                        px[4] * x31 + (1.0 / 2.0) * px[5] * x18 * x3 + (1.0 / 2.0) * px[6] * x18 * x8 - px[7] * x32;
                T x38 = -py[0] * x5 - py[1] * x7 - py[2] * x10 - py[3] * x11 + (1.0 / 2.0) * py[4] * x1 * x3 +
                        (1.0 / 2.0) * py[5] * x3 * x6 + (1.0 / 2.0) * py[6] * x6 * x8 + (1.0 / 2.0) * py[7] * x1 * x8;
                T x39 = x37 * x38;
                T x40 = x12 * x22 * x35 - x22 * x39 + x23 * x33 * x38 - x23 * x36 + x24 * x28 * x37 - x28 * x34;
                T x41 = 1.0 / x40;
                T x42 = x14 * x41;
                T x43 = x27 * x42;
                T x44 = -x22 * x38 + x24 * x28;
                T x45 = (1.0 / 4.0) * y;
                T x46 = x45 - 0.25;
                T x47 = x42 * x46;
                T x48 = -x12 * x28 + x23 * x38;
                T x49 = x41 * x5;
                T x50 = x25 * x43 + x44 * x47 - x48 * x49;
                T x51 = x24 * x37 - x34;
                T x52 = x33 * x38 - x36;
                T x53 = x12 * x35 - x39;
                T x54 = x43 * x51 + x47 * x52 - x49 * x53;
                T x55 = -x22 * x37 + x23 * x33;
                T x56 = x22 * x35 - x28 * x33;
                T x57 = -x23 * x35 + x28 * x37;
                T x58 = x43 * x55 + x47 * x56 - x49 * x57;
                T x59 = weight * x40;
                T x60 = -x42 * x46;
                T x61 = x26 + 0.25;
                T x62 = -x42 * x61;
                T x63 = x41 * x7;
                T x64 = x25 * x62 + x44 * x60 - x48 * x63;
                T x65 = x51 * x62 + x52 * x60 - x53 * x63;
                T x66 = x55 * x62 + x56 * x60 - x57 * x63;
                T x67 = x59 * (x50 * x64 + x54 * x65 + x58 * x66);
                T x68 = x42 * x61;
                T x69 = x45 + 0.25;
                T x70 = x42 * x69;
                T x71 = x10 * x41;
                T x72 = x48 * x71;
                T x73 = x25 * x68 + x44 * x70 - x72;
                T x74 = x53 * x71;
                T x75 = x51 * x68 + x52 * x70 - x74;
                T x76 = x57 * x71;
                T x77 = x55 * x68 + x56 * x70 - x76;
                T x78 = x59 * (x50 * x73 + x54 * x75 + x58 * x77);
                T x79 = -x27;
                T x80 = x42 * x79;
                T x81 = -x42 * x69;
                T x82 = x11 * x41;
                T x83 = x25 * x80 + x44 * x81 - x48 * x82;
                T x84 = x51 * x80 + x52 * x81 - x53 * x82;
                T x85 = x55 * x80 + x56 * x81 - x57 * x82;
                T x86 = x59 * (x50 * x83 + x54 * x84 + x58 * x85);
                T x87 = x3 * x41;
                T x88 = x79 * x87;
                T x89 = (1.0 / 4.0) * z + 0.25;
                T x90 = -x87 * x89;
                T x91 = x20 * x41;
                T x92 = -x25 * x91 + x44 * x90 + x48 * x88;
                T x93 = -x51 * x91 + x52 * x90 + x53 * x88;
                T x94 = -x55 * x91 + x56 * x90 + x57 * x88;
                T x95 = x59 * (x50 * x92 + x54 * x93 + x58 * x94);
                T x96 = x61 * x87;
                T x97 = x87 * x89;
                T x98 = x21 * x41;
                T x99 = -x25 * x98 + x44 * x97 + x48 * x96;
                T x100 = -x51 * x98 + x52 * x97 + x53 * x96;
                T x101 = -x55 * x98 + x56 * x97 + x57 * x96;
                T x102 = x59 * (x100 * x54 + x101 * x58 + x50 * x99);
                T x103 = x18 * x41;
                T x104 = x103 * x61;
                T x105 = x103 * x69;
                T x106 = x104 * x25 + x105 * x44 + x72;
                T x107 = x104 * x51 + x105 * x52 + x74;
                T x108 = x104 * x55 + x105 * x56 + x76;
                T x109 = x59 * (x106 * x50 + x107 * x54 + x108 * x58);
                T x110 = x1 * x41;
                T x111 = x110 * x69;
                T x112 = x110 * x89;
                T x113 = x32 * x41;
                T x114 = x111 * x48 + x112 * x25 - x113 * x44;
                T x115 = x111 * x53 + x112 * x51 - x113 * x52;
                T x116 = x111 * x57 + x112 * x55 - x113 * x56;
                T x117 = x59 * (x114 * x50 + x115 * x54 + x116 * x58);
                T x118 = x59 * (x64 * x73 + x65 * x75 + x66 * x77);
                T x119 = x59 * (x64 * x83 + x65 * x84 + x66 * x85);
                T x120 = x59 * (x64 * x92 + x65 * x93 + x66 * x94);
                T x121 = x59 * (x100 * x65 + x101 * x66 + x64 * x99);
                T x122 = x59 * (x106 * x64 + x107 * x65 + x108 * x66);
                T x123 = x59 * (x114 * x64 + x115 * x65 + x116 * x66);
                T x124 = x59 * (x73 * x83 + x75 * x84 + x77 * x85);
                T x125 = x59 * (x73 * x92 + x75 * x93 + x77 * x94);
                T x126 = x59 * (x100 * x75 + x101 * x77 + x73 * x99);
                T x127 = x59 * (x106 * x73 + x107 * x75 + x108 * x77);
                T x128 = x59 * (x114 * x73 + x115 * x75 + x116 * x77);
                T x129 = x59 * (x83 * x92 + x84 * x93 + x85 * x94);
                T x130 = x59 * (x100 * x84 + x101 * x85 + x83 * x99);
                T x131 = x59 * (x106 * x83 + x107 * x84 + x108 * x85);
                T x132 = x59 * (x114 * x83 + x115 * x84 + x116 * x85);
                T x133 = x59 * (x100 * x93 + x101 * x94 + x92 * x99);
                T x134 = x59 * (x106 * x92 + x107 * x93 + x108 * x94);
                T x135 = x59 * (x114 * x92 + x115 * x93 + x116 * x94);
                T x136 = x59 * (x100 * x107 + x101 * x108 + x106 * x99);
                T x137 = x59 * (x100 * x115 + x101 * x116 + x114 * x99);
                T x138 = x59 * (x106 * x114 + x107 * x115 + x108 * x116);
                T x139 = u[0] * x50 + u[1] * x64 + u[2] * x73 + u[3] * x83 + u[4] * x92 + u[5] * x99 + u[6] * x106 +
                         u[7] * x114;
                T x140 = u[0] * x54 + u[1] * x65 + u[2] * x75 + u[3] * x84 + u[4] * x93 + u[5] * x100 + u[6] * x107 +
                         u[7] * x115;
                T x141 = u[0] * x58 + u[1] * x66 + u[2] * x77 + u[3] * x85 + u[4] * x94 + u[5] * x101 + u[6] * x108 +
                         u[7] * x116;
                T x142 = 8 * x59;
                H[0] += x59 * (pow(x50, 2) + pow(x54, 2) + pow(x58, 2));
                H[1] += x67;
                H[2] += x78;
                H[3] += x86;
                H[4] += x95;
                H[5] += x102;
                H[6] += x109;
                H[7] += x117;
                H[8] += x67;
                H[9] += x59 * (pow(x64, 2) + pow(x65, 2) + pow(x66, 2));
                H[10] += x118;
                H[11] += x119;
                H[12] += x120;
                H[13] += x121;
                H[14] += x122;
                H[15] += x123;
                H[16] += x78;
                H[17] += x118;
                H[18] += x59 * (pow(x73, 2) + pow(x75, 2) + pow(x77, 2));
                H[19] += x124;
                H[20] += x125;
                H[21] += x126;
                H[22] += x127;
                H[23] += x128;
                H[24] += x86;
                H[25] += x119;
                H[26] += x124;
                H[27] += x59 * (pow(x83, 2) + pow(x84, 2) + pow(x85, 2));
                H[28] += x129;
                H[29] += x130;
                H[30] += x131;
                H[31] += x132;
                H[32] += x95;
                H[33] += x120;
                H[34] += x125;
                H[35] += x129;
                H[36] += x59 * (pow(x92, 2) + pow(x93, 2) + pow(x94, 2));
                H[37] += x133;
                H[38] += x134;
                H[39] += x135;
                H[40] += x102;
                H[41] += x121;
                H[42] += x126;
                H[43] += x130;
                H[44] += x133;
                H[45] += x59 * (pow(x100, 2) + pow(x101, 2) + pow(x99, 2));
                H[46] += x136;
                H[47] += x137;
                H[48] += x109;
                H[49] += x122;
                H[50] += x127;
                H[51] += x131;
                H[52] += x134;
                H[53] += x136;
                H[54] += x59 * (pow(x106, 2) + pow(x107, 2) + pow(x108, 2));
                H[55] += x138;
                H[56] += x117;
                H[57] += x123;
                H[58] += x128;
                H[59] += x132;
                H[60] += x135;
                H[61] += x137;
                H[62] += x138;
                H[63] += x59 * (pow(x114, 2) + pow(x115, 2) + pow(x116, 2));
                g[0] += x142 * (x139 * x50 + x140 * x54 + x141 * x58);
                g[1] += x142 * (x139 * x64 + x140 * x65 + x141 * x66);
                g[2] += x142 * (x139 * x73 + x140 * x75 + x141 * x77);
                g[3] += x142 * (x139 * x83 + x140 * x84 + x141 * x85);
                g[4] += x142 * (x139 * x92 + x140 * x93 + x141 * x94);
                g[5] += x142 * (x100 * x140 + x101 * x141 + x139 * x99);
                g[6] += x142 * (x106 * x139 + x107 * x140 + x108 * x141);
                g[7] += x142 * (x114 * x139 + x115 * x140 + x116 * x141);
                e += x59 * (pow(x139, 2) + pow(x140, 2) + pow(x141, 2));
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FE>
        using LaplaceOperatorHex8 = utopia::kokkos::AutoKernel<
            FE,
            utopia::kernels::LaplaceOperator<utopia::kernels::Hex8<typename FE::Scalar, typename FE::Scalar>>,
            3>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_LaplaceOperator_Hex8_3_IMPL_hpp
