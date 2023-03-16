#ifndef UTOPIA_TPL_MATERIAL_LaplaceOperator_Pentatope5_4_IMPL_hpp
#define UTOPIA_TPL_MATERIAL_LaplaceOperator_Pentatope5_4_IMPL_hpp

#include "utopia_Algorithms.hpp"
#include "utopia_Input.hpp"

#include "utopia_kokkos_AutoKernel.hpp"

#include "utopia_fe_Pentatope5_4.hpp"
#include "utopia_material_LaplaceOperator.hpp"

#ifndef UTOPIA_RESTRICT
#define UTOPIA_RESTRICT __restrict__
#endif

namespace utopia {
    namespace kernels {

        /**
         * Specialization of LaplaceOperator for symmetric element pair trial=test=Pentatope5
         */
        template <typename T, typename GeoT>
        class LaplaceOperator<Pentatope5<T, GeoT>> {
        public:
            using ElemT = Pentatope5<T, GeoT>;
            static constexpr int Dim = ElemT::Dim;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "LaplaceOperator<Pentatope5>"; }

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
                const GeoT *UTOPIA_RESTRICT pt,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T t,
                const T weight,
                T *UTOPIA_RESTRICT H) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 5*ADD + 25*ADDAUGMENTEDASSIGNMENT + 20*MUL + 19*POW
                //	- Subexpressions: 67*ADD + DIV + 203*MUL + 8*NEG + POW + 64*SUB
                T x0 = -pt[0] + pt[1];
                T x1 = -py[0] + py[3];
                T x2 = -pz[0] + pz[2];
                T x3 = x1 * x2;
                T x4 = py[0] - py[1];
                T x5 = -x4;
                T x6 = -pz[0] + pz[3];
                T x7 = -pt[0] + pt[2];
                T x8 = x6 * x7;
                T x9 = -py[0] + py[2];
                T x10 = -pz[0] + pz[1];
                T x11 = -pt[0] + pt[3];
                T x12 = x10 * x11;
                T x13 = x11 * x2;
                T x14 = x6 * x9;
                T x15 = x1 * x7;
                T x16 = -x0 * x14 + x0 * x3 - x10 * x15 + x12 * x9 + x13 * x4 + x5 * x8;
                T x17 = px[0] - px[2];
                T x18 = -x17;
                T x19 = -py[0] + py[4];
                T x20 = x19 * x6;
                T x21 = -px[0] + px[3];
                T x22 = -pz[0] + pz[4];
                T x23 = x22 * x9;
                T x24 = -px[0] + px[4];
                T x25 = x1 * x22;
                T x26 = x19 * x2;
                T x27 = -x14 * x24 + x17 * x25 + x18 * x20 + x21 * x23 - x21 * x26 + x24 * x3;
                T x28 = x24 * x7;
                T x29 = x11 * x18;
                T x30 = -pt[0] + pt[4];
                T x31 = x21 * x30;
                T x32 = x30 * x6;
                T x33 = x21 * x7;
                T x34 = x11 * x24;
                T x35 = x17 * x32 + x2 * x31 - x2 * x34 + x22 * x29 - x22 * x33 + x28 * x6;
                T x36 = px[0] - px[1];
                T x37 = -x36;
                T x38 = -x11 * x23 + x11 * x26 + x14 * x30 - x20 * x7 + x25 * x7 - x3 * x30;
                T x39 = x1 * x30;
                T x40 = -x1 * x28 + x18 * x39 - x19 * x29 + x19 * x33 - x31 * x9 + x34 * x9;
                T x41 = x0 * x27 + x10 * x40 + x35 * x5 + x37 * x38;
                T x42 = 1.0 / x41;
                T x43 = x16 * x42;
                T x44 = x11 * x22;
                T x45 = x0 * x20 - x0 * x25 + x10 * x39 - x12 * x19 + x32 * x4 + x44 * x5;
                T x46 = x42 * x45;
                T x47 = x38 * x42;
                T x48 = -x43 - x46 - x47;
                T x49 = x0 * x9;
                T x50 = x0 * x18;
                T x51 = -x1 * x50 + x11 * x36 * x9 + x15 * x37 + x21 * x49 + x29 * x5 - x33 * x5;
                T x52 = x42 * x51;
                T x53 = x0 * x24;
                T x54 = x19 * x37;
                T x55 = x0 * x21;
                T x56 = x1 * x53 + x11 * x54 - x19 * x55 + x31 * x5 - x34 * x5 + x36 * x39;
                T x57 = x42 * x56;
                T x58 = x18 * x30;
                T x59 = x19 * x50 - x24 * x49 + x28 * x5 + x30 * x37 * x9 - x5 * x58 - x54 * x7;
                T x60 = x42 * x59;
                T x61 = x40 * x42;
                T x62 = -x52 - x57 - x60 - x61;
                T x63 = -x10 * x28 + x10 * x58 + x2 * x30 * x36 + x2 * x53 + x22 * x37 * x7 - x22 * x50;
                T x64 = x42 * x63;
                T x65 = x35 * x42;
                T x66 = -x10 * x29 + x10 * x33 + x13 * x37 - x2 * x55 - x37 * x8 + x50 * x6;
                T x67 = x42 * x66;
                T x68 = -x10 * x31 + x10 * x34 + x22 * x55 + x32 * x37 - x37 * x44 - x53 * x6;
                T x69 = x42 * x68;
                T x70 = -x64 - x65 - x67 - x69;
                T x71 = x18 * x5;
                T x72 = x10 * x9;
                T x73 = x10 * x18;
                T x74 = x2 * x5;
                T x75 = -x19 * x73 + x22 * x71 + x23 * x36 + x24 * x72 - x24 * x74 + x26 * x37;
                T x76 = x42 * x75;
                T x77 = x27 * x42;
                T x78 = x1 * x73 + x14 * x37 - x21 * x72 + x21 * x74 - x3 * x37 - x6 * x71;
                T x79 = x42 * x78;
                T x80 = -x1 * x10 * x24 + x10 * x19 * x21 - x20 * x37 - x21 * x22 * x5 + x24 * x5 * x6 + x25 * x37;
                T x81 = x42 * x80;
                T x82 = -x76 - x77 - x79 - x81;
                T x83 = weight * x41;
                T x84 = x83 * (x47 * x48 + x61 * x62 + x65 * x70 + x77 * x82);
                T x85 = x83 * (x46 * x48 + x57 * x62 + x69 * x70 + x81 * x82);
                T x86 = x83 * (x60 * x62 + x64 * x70 + x76 * x82);
                T x87 = x83 * (x43 * x48 + x52 * x62 + x67 * x70 + x79 * x82);
                T x88 = pow(x41, -2);
                T x89 = x40 * x88;
                T x90 = x35 * x88;
                T x91 = x38 * x88;
                T x92 = x27 * x88;
                T x93 = x83 * (x45 * x91 + x56 * x89 + x68 * x90 + x80 * x92);
                T x94 = x83 * (x59 * x89 + x63 * x90 + x75 * x92);
                T x95 = x83 * (x16 * x91 + x51 * x89 + x66 * x90 + x78 * x92);
                T x96 = x56 * x88;
                T x97 = x68 * x88;
                T x98 = x80 * x88;
                T x99 = x83 * (x59 * x96 + x63 * x97 + x75 * x98);
                T x100 = x83 * (x16 * x45 * x88 + x51 * x96 + x66 * x97 + x78 * x98);
                T x101 = x83 * (x51 * x59 * x88 + x63 * x66 * x88 + x75 * x78 * x88);
                H[0] += x83 * (pow(x48, 2) + pow(x62, 2) + pow(x70, 2) + pow(x82, 2));
                H[1] += x84;
                H[2] += x85;
                H[3] += x86;
                H[4] += x87;
                H[5] += x84;
                H[6] += x83 * (pow(x27, 2) * x88 + pow(x35, 2) * x88 + pow(x38, 2) * x88 + pow(x40, 2) * x88);
                H[7] += x93;
                H[8] += x94;
                H[9] += x95;
                H[10] += x85;
                H[11] += x93;
                H[12] += x83 * (pow(x45, 2) * x88 + pow(x56, 2) * x88 + pow(x68, 2) * x88 + pow(x80, 2) * x88);
                H[13] += x99;
                H[14] += x100;
                H[15] += x86;
                H[16] += x94;
                H[17] += x99;
                H[18] += x83 * (pow(x59, 2) * x88 + pow(x63, 2) * x88 + pow(x75, 2) * x88);
                H[19] += x101;
                H[20] += x87;
                H[21] += x95;
                H[22] += x100;
                H[23] += x101;
                H[24] += x83 * (pow(x16, 2) * x88 + pow(x51, 2) * x88 + pow(x66, 2) * x88 + pow(x78, 2) * x88);
            }

            UTOPIA_FUNCTION void apply(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                const GeoT *UTOPIA_RESTRICT pt,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T t,
                const T weight,
                T *UTOPIA_RESTRICT Hx) const {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                UTOPIA_UNUSED(z);
                UTOPIA_UNUSED(t);
                // FLOATING POINT OPS!
                //	- Result: 5*ADD + 5*ADDAUGMENTEDASSIGNMENT + 24*MUL
                //	- Subexpressions: 56*ADD + DIV + 166*MUL + 7*NEG + 64*SUB
                T x0 = -pt[0] + pt[1];
                T x1 = -py[0] + py[3];
                T x2 = -pz[0] + pz[2];
                T x3 = x1 * x2;
                T x4 = py[0] - py[1];
                T x5 = -x4;
                T x6 = -pz[0] + pz[3];
                T x7 = -pt[0] + pt[2];
                T x8 = x6 * x7;
                T x9 = -py[0] + py[2];
                T x10 = -pz[0] + pz[1];
                T x11 = -pt[0] + pt[3];
                T x12 = x10 * x11;
                T x13 = x11 * x2;
                T x14 = x6 * x9;
                T x15 = x1 * x7;
                T x16 = px[0] - px[2];
                T x17 = -x16;
                T x18 = -py[0] + py[4];
                T x19 = x18 * x6;
                T x20 = -px[0] + px[3];
                T x21 = -pz[0] + pz[4];
                T x22 = x21 * x9;
                T x23 = -px[0] + px[4];
                T x24 = x1 * x21;
                T x25 = x18 * x2;
                T x26 = -x14 * x23 + x16 * x24 + x17 * x19 + x20 * x22 - x20 * x25 + x23 * x3;
                T x27 = x23 * x7;
                T x28 = x11 * x17;
                T x29 = -pt[0] + pt[4];
                T x30 = x20 * x29;
                T x31 = x29 * x6;
                T x32 = x20 * x7;
                T x33 = x11 * x23;
                T x34 = x16 * x31 + x2 * x30 - x2 * x33 + x21 * x28 - x21 * x32 + x27 * x6;
                T x35 = px[0] - px[1];
                T x36 = -x35;
                T x37 = -x11 * x22 + x11 * x25 + x14 * x29 - x19 * x7 + x24 * x7 - x29 * x3;
                T x38 = x1 * x29;
                T x39 = -x1 * x27 + x17 * x38 - x18 * x28 + x18 * x32 - x30 * x9 + x33 * x9;
                T x40 = x0 * x26 + x10 * x39 + x34 * x5 + x36 * x37;
                T x41 = 1.0 / x40;
                T x42 = x41 * (-x0 * x14 + x0 * x3 - x10 * x15 + x12 * x9 + x13 * x4 + x5 * x8);
                T x43 = x11 * x21;
                T x44 = x41 * (x0 * x19 - x0 * x24 + x10 * x38 - x12 * x18 + x31 * x4 + x43 * x5);
                T x45 = x37 * x41;
                T x46 = -x42 - x44 - x45;
                T x47 = u[0] * x46 + u[1] * x45 + u[2] * x44 + u[4] * x42;
                T x48 = x0 * x9;
                T x49 = x0 * x17;
                T x50 = x41 * (-x1 * x49 + x11 * x35 * x9 + x15 * x36 + x20 * x48 + x28 * x5 - x32 * x5);
                T x51 = x0 * x23;
                T x52 = x18 * x36;
                T x53 = x0 * x20;
                T x54 = x41 * (x1 * x51 + x11 * x52 - x18 * x53 + x30 * x5 - x33 * x5 + x35 * x38);
                T x55 = x17 * x29;
                T x56 = x41 * (x18 * x49 - x23 * x48 + x27 * x5 + x29 * x36 * x9 - x5 * x55 - x52 * x7);
                T x57 = x39 * x41;
                T x58 = -x50 - x54 - x56 - x57;
                T x59 = u[0] * x58 + u[1] * x57 + u[2] * x54 + u[3] * x56 + u[4] * x50;
                T x60 = x41 * (-x10 * x27 + x10 * x55 + x2 * x29 * x35 + x2 * x51 + x21 * x36 * x7 - x21 * x49);
                T x61 = x34 * x41;
                T x62 = x41 * (-x10 * x28 + x10 * x32 + x13 * x36 - x2 * x53 - x36 * x8 + x49 * x6);
                T x63 = x41 * (-x10 * x30 + x10 * x33 + x21 * x53 + x31 * x36 - x36 * x43 - x51 * x6);
                T x64 = -x60 - x61 - x62 - x63;
                T x65 = u[0] * x64 + u[1] * x61 + u[2] * x63 + u[3] * x60 + u[4] * x62;
                T x66 = x17 * x5;
                T x67 = x10 * x9;
                T x68 = x10 * x17;
                T x69 = x2 * x5;
                T x70 = x41 * (-x18 * x68 + x21 * x66 + x22 * x35 + x23 * x67 - x23 * x69 + x25 * x36);
                T x71 = x26 * x41;
                T x72 = x41 * (x1 * x68 + x14 * x36 - x20 * x67 + x20 * x69 - x3 * x36 - x6 * x66);
                T x73 =
                    x41 * (-x1 * x10 * x23 + x10 * x18 * x20 - x19 * x36 - x20 * x21 * x5 + x23 * x5 * x6 + x24 * x36);
                T x74 = -x70 - x71 - x72 - x73;
                T x75 = u[0] * x74 + u[1] * x71 + u[2] * x73 + u[3] * x70 + u[4] * x72;
                T x76 = 5 * weight * x40;
                Hx[0] += x76 * (x46 * x47 + x58 * x59 + x64 * x65 + x74 * x75);
                Hx[1] += x76 * (x45 * x47 + x57 * x59 + x61 * x65 + x71 * x75);
                Hx[2] += x76 * (x44 * x47 + x54 * x59 + x63 * x65 + x73 * x75);
                Hx[3] += x76 * (x56 * x59 + x60 * x65 + x70 * x75);
                Hx[4] += x76 * (x42 * x47 + x50 * x59 + x62 * x65 + x72 * x75);
            }

            UTOPIA_FUNCTION void gradient(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                const GeoT *UTOPIA_RESTRICT pt,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T t,
                const T weight,
                T *UTOPIA_RESTRICT g) const {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                UTOPIA_UNUSED(z);
                UTOPIA_UNUSED(t);
                // FLOATING POINT OPS!
                //	- Result: 5*ADD + 5*ADDAUGMENTEDASSIGNMENT + 24*MUL
                //	- Subexpressions: 56*ADD + DIV + 166*MUL + 7*NEG + 64*SUB
                T x0 = -pt[0] + pt[1];
                T x1 = -py[0] + py[3];
                T x2 = -pz[0] + pz[2];
                T x3 = x1 * x2;
                T x4 = py[0] - py[1];
                T x5 = -x4;
                T x6 = -pz[0] + pz[3];
                T x7 = -pt[0] + pt[2];
                T x8 = x6 * x7;
                T x9 = -py[0] + py[2];
                T x10 = -pz[0] + pz[1];
                T x11 = -pt[0] + pt[3];
                T x12 = x10 * x11;
                T x13 = x11 * x2;
                T x14 = x6 * x9;
                T x15 = x1 * x7;
                T x16 = px[0] - px[2];
                T x17 = -x16;
                T x18 = -py[0] + py[4];
                T x19 = x18 * x6;
                T x20 = -px[0] + px[3];
                T x21 = -pz[0] + pz[4];
                T x22 = x21 * x9;
                T x23 = -px[0] + px[4];
                T x24 = x1 * x21;
                T x25 = x18 * x2;
                T x26 = -x14 * x23 + x16 * x24 + x17 * x19 + x20 * x22 - x20 * x25 + x23 * x3;
                T x27 = x23 * x7;
                T x28 = x11 * x17;
                T x29 = -pt[0] + pt[4];
                T x30 = x20 * x29;
                T x31 = x29 * x6;
                T x32 = x20 * x7;
                T x33 = x11 * x23;
                T x34 = x16 * x31 + x2 * x30 - x2 * x33 + x21 * x28 - x21 * x32 + x27 * x6;
                T x35 = px[0] - px[1];
                T x36 = -x35;
                T x37 = -x11 * x22 + x11 * x25 + x14 * x29 - x19 * x7 + x24 * x7 - x29 * x3;
                T x38 = x1 * x29;
                T x39 = -x1 * x27 + x17 * x38 - x18 * x28 + x18 * x32 - x30 * x9 + x33 * x9;
                T x40 = x0 * x26 + x10 * x39 + x34 * x5 + x36 * x37;
                T x41 = 1.0 / x40;
                T x42 = x41 * (-x0 * x14 + x0 * x3 - x10 * x15 + x12 * x9 + x13 * x4 + x5 * x8);
                T x43 = x11 * x21;
                T x44 = x41 * (x0 * x19 - x0 * x24 + x10 * x38 - x12 * x18 + x31 * x4 + x43 * x5);
                T x45 = x37 * x41;
                T x46 = -x42 - x44 - x45;
                T x47 = u[0] * x46 + u[1] * x45 + u[2] * x44 + u[4] * x42;
                T x48 = x0 * x9;
                T x49 = x0 * x17;
                T x50 = x41 * (-x1 * x49 + x11 * x35 * x9 + x15 * x36 + x20 * x48 + x28 * x5 - x32 * x5);
                T x51 = x0 * x23;
                T x52 = x18 * x36;
                T x53 = x0 * x20;
                T x54 = x41 * (x1 * x51 + x11 * x52 - x18 * x53 + x30 * x5 - x33 * x5 + x35 * x38);
                T x55 = x17 * x29;
                T x56 = x41 * (x18 * x49 - x23 * x48 + x27 * x5 + x29 * x36 * x9 - x5 * x55 - x52 * x7);
                T x57 = x39 * x41;
                T x58 = -x50 - x54 - x56 - x57;
                T x59 = u[0] * x58 + u[1] * x57 + u[2] * x54 + u[3] * x56 + u[4] * x50;
                T x60 = x41 * (-x10 * x27 + x10 * x55 + x2 * x29 * x35 + x2 * x51 + x21 * x36 * x7 - x21 * x49);
                T x61 = x34 * x41;
                T x62 = x41 * (-x10 * x28 + x10 * x32 + x13 * x36 - x2 * x53 - x36 * x8 + x49 * x6);
                T x63 = x41 * (-x10 * x30 + x10 * x33 + x21 * x53 + x31 * x36 - x36 * x43 - x51 * x6);
                T x64 = -x60 - x61 - x62 - x63;
                T x65 = u[0] * x64 + u[1] * x61 + u[2] * x63 + u[3] * x60 + u[4] * x62;
                T x66 = x17 * x5;
                T x67 = x10 * x9;
                T x68 = x10 * x17;
                T x69 = x2 * x5;
                T x70 = x41 * (-x18 * x68 + x21 * x66 + x22 * x35 + x23 * x67 - x23 * x69 + x25 * x36);
                T x71 = x26 * x41;
                T x72 = x41 * (x1 * x68 + x14 * x36 - x20 * x67 + x20 * x69 - x3 * x36 - x6 * x66);
                T x73 =
                    x41 * (-x1 * x10 * x23 + x10 * x18 * x20 - x19 * x36 - x20 * x21 * x5 + x23 * x5 * x6 + x24 * x36);
                T x74 = -x70 - x71 - x72 - x73;
                T x75 = u[0] * x74 + u[1] * x71 + u[2] * x73 + u[3] * x70 + u[4] * x72;
                T x76 = 5 * weight * x40;
                g[0] += x76 * (x46 * x47 + x58 * x59 + x64 * x65 + x74 * x75);
                g[1] += x76 * (x45 * x47 + x57 * x59 + x61 * x65 + x71 * x75);
                g[2] += x76 * (x44 * x47 + x54 * x59 + x63 * x65 + x73 * x75);
                g[3] += x76 * (x56 * x59 + x60 * x65 + x70 * x75);
                g[4] += x76 * (x42 * x47 + x50 * x59 + x62 * x65 + x72 * x75);
            }

            UTOPIA_FUNCTION void value(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                const GeoT *UTOPIA_RESTRICT pt,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T t,
                const T weight,
                T &e) const {
                using namespace utopia::device;
                // Automatically generated
                // FLOATING POINT OPS!
                //	- Result: 9*ADD + ADDAUGMENTEDASSIGNMENT + 35*MUL + 4*POW
                //	- Subexpressions: 41*ADD + DIV + 145*MUL + 3*NEG + 53*SUB
                T x0 = -pt[0] + pt[1];
                T x1 = px[0] - px[2];
                T x2 = -x1;
                T x3 = -py[0] + py[4];
                T x4 = -pz[0] + pz[3];
                T x5 = x3 * x4;
                T x6 = -px[0] + px[3];
                T x7 = -pz[0] + pz[4];
                T x8 = -py[0] + py[2];
                T x9 = x7 * x8;
                T x10 = -px[0] + px[4];
                T x11 = -py[0] + py[3];
                T x12 = -pz[0] + pz[2];
                T x13 = x11 * x12;
                T x14 = x11 * x7;
                T x15 = x12 * x3;
                T x16 = x4 * x8;
                T x17 = x1 * x14 + x10 * x13 - x10 * x16 - x15 * x6 + x2 * x5 + x6 * x9;
                T x18 = py[0] - py[1];
                T x19 = -x18;
                T x20 = -pt[0] + pt[2];
                T x21 = x10 * x20;
                T x22 = -pt[0] + pt[3];
                T x23 = x2 * x22;
                T x24 = -pt[0] + pt[4];
                T x25 = x24 * x6;
                T x26 = x24 * x4;
                T x27 = x20 * x6;
                T x28 = x10 * x22;
                T x29 = x1 * x26 + x12 * x25 - x12 * x28 + x21 * x4 + x23 * x7 - x27 * x7;
                T x30 = px[0] - px[1];
                T x31 = -x30;
                T x32 = -x13 * x24 + x14 * x20 + x15 * x22 + x16 * x24 - x20 * x5 - x22 * x9;
                T x33 = -pz[0] + pz[1];
                T x34 = x11 * x24;
                T x35 = -x11 * x21 + x2 * x34 - x23 * x3 - x25 * x8 + x27 * x3 + x28 * x8;
                T x36 = x0 * x17 + x19 * x29 + x31 * x32 + x33 * x35;
                T x37 = x22 * x7;
                T x38 = x22 * x33;
                T x39 = 1.0 / x36;
                T x40 = x39 * (-x0 * x14 + x0 * x5 + x18 * x26 + x19 * x37 - x3 * x38 + x33 * x34);
                T x41 = x20 * x4;
                T x42 = x12 * x22;
                T x43 = x11 * x20;
                T x44 = x39 * (x0 * x13 - x0 * x16 + x18 * x42 + x19 * x41 - x33 * x43 + x38 * x8);
                T x45 = x32 * x39;
                T x46 = x0 * x10;
                T x47 = x3 * x31;
                T x48 = x0 * x6;
                T x49 = x39 * (x11 * x46 + x19 * x25 - x19 * x28 + x22 * x47 - x3 * x48 + x30 * x34);
                T x50 = x0 * x8;
                T x51 = x0 * x2;
                T x52 = x39 * (-x11 * x51 + x19 * x23 - x19 * x27 + x22 * x30 * x8 + x31 * x43 + x50 * x6);
                T x53 = x35 * x39;
                T x54 = x2 * x24;
                T x55 = x39 * (-x10 * x50 + x19 * x21 - x19 * x54 - x20 * x47 + x24 * x31 * x8 + x3 * x51);
                T x56 = x29 * x39;
                T x57 = x39 * (x12 * x24 * x30 + x12 * x46 + x20 * x31 * x7 - x21 * x33 + x33 * x54 - x51 * x7);
                T x58 = x39 * (-x25 * x33 + x26 * x31 + x28 * x33 - x31 * x37 - x4 * x46 + x48 * x7);
                T x59 = x39 * (-x12 * x48 - x23 * x33 + x27 * x33 - x31 * x41 + x31 * x42 + x4 * x51);
                T x60 = x17 * x39;
                T x61 = x19 * x2;
                T x62 = x33 * x8;
                T x63 = x2 * x33;
                T x64 = x12 * x19;
                T x65 = x39 * (x10 * x62 - x10 * x64 + x15 * x31 - x3 * x63 + x30 * x9 + x61 * x7);
                T x66 =
                    x39 * (-x10 * x11 * x33 + x10 * x19 * x4 + x14 * x31 - x19 * x6 * x7 + x3 * x33 * x6 - x31 * x5);
                T x67 = x39 * (x11 * x63 - x13 * x31 + x16 * x31 - x4 * x61 - x6 * x62 + x6 * x64);
                e += weight * x36 *
                     (pow(u[0] * (-x40 - x44 - x45) + u[1] * x45 + u[2] * x40 + u[4] * x44, 2) +
                      pow(u[0] * (-x49 - x52 - x53 - x55) + u[1] * x53 + u[2] * x49 + u[3] * x55 + u[4] * x52, 2) +
                      pow(u[0] * (-x56 - x57 - x58 - x59) + u[1] * x56 + u[2] * x58 + u[3] * x57 + u[4] * x59, 2) +
                      pow(u[0] * (-x60 - x65 - x66 - x67) + u[1] * x60 + u[2] * x66 + u[3] * x65 + u[4] * x67, 2));
            }

            UTOPIA_FUNCTION void eval(
                // Element coordinates
                const GeoT *UTOPIA_RESTRICT px,
                const GeoT *UTOPIA_RESTRICT py,
                const GeoT *UTOPIA_RESTRICT pz,
                const GeoT *UTOPIA_RESTRICT pt,
                // Coefficients
                const T *UTOPIA_RESTRICT u,
                // Quadrature rule
                const T x,
                const T y,
                const T z,
                const T t,
                const T weight,
                T &e,
                T *UTOPIA_RESTRICT H,
                T *UTOPIA_RESTRICT g) const {
                using namespace utopia::device;
                // Automatically generated

                // Unused variables
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                UTOPIA_UNUSED(z);
                UTOPIA_UNUSED(t);
                // FLOATING POINT OPS!
                //	- Result: 11*ADD + 31*ADDAUGMENTEDASSIGNMENT + 45*MUL + 23*POW
                //	- Subexpressions: 82*ADD + DIV + 223*MUL + 8*NEG + POW + 64*SUB
                T x0 = -pt[0] + pt[1];
                T x1 = -py[0] + py[3];
                T x2 = -pz[0] + pz[2];
                T x3 = x1 * x2;
                T x4 = py[0] - py[1];
                T x5 = -x4;
                T x6 = -pz[0] + pz[3];
                T x7 = -pt[0] + pt[2];
                T x8 = x6 * x7;
                T x9 = -py[0] + py[2];
                T x10 = -pz[0] + pz[1];
                T x11 = -pt[0] + pt[3];
                T x12 = x10 * x11;
                T x13 = x11 * x2;
                T x14 = x6 * x9;
                T x15 = x1 * x7;
                T x16 = -x0 * x14 + x0 * x3 - x10 * x15 + x12 * x9 + x13 * x4 + x5 * x8;
                T x17 = px[0] - px[2];
                T x18 = -x17;
                T x19 = -py[0] + py[4];
                T x20 = x19 * x6;
                T x21 = -px[0] + px[3];
                T x22 = -pz[0] + pz[4];
                T x23 = x22 * x9;
                T x24 = -px[0] + px[4];
                T x25 = x1 * x22;
                T x26 = x19 * x2;
                T x27 = -x14 * x24 + x17 * x25 + x18 * x20 + x21 * x23 - x21 * x26 + x24 * x3;
                T x28 = x24 * x7;
                T x29 = x11 * x18;
                T x30 = -pt[0] + pt[4];
                T x31 = x21 * x30;
                T x32 = x30 * x6;
                T x33 = x21 * x7;
                T x34 = x11 * x24;
                T x35 = x17 * x32 + x2 * x31 - x2 * x34 + x22 * x29 - x22 * x33 + x28 * x6;
                T x36 = px[0] - px[1];
                T x37 = -x36;
                T x38 = -x11 * x23 + x11 * x26 + x14 * x30 - x20 * x7 + x25 * x7 - x3 * x30;
                T x39 = x1 * x30;
                T x40 = -x1 * x28 + x18 * x39 - x19 * x29 + x19 * x33 - x31 * x9 + x34 * x9;
                T x41 = x0 * x27 + x10 * x40 + x35 * x5 + x37 * x38;
                T x42 = 1.0 / x41;
                T x43 = x16 * x42;
                T x44 = x11 * x22;
                T x45 = x0 * x20 - x0 * x25 + x10 * x39 - x12 * x19 + x32 * x4 + x44 * x5;
                T x46 = x42 * x45;
                T x47 = x38 * x42;
                T x48 = -x43 - x46 - x47;
                T x49 = x0 * x9;
                T x50 = x0 * x18;
                T x51 = -x1 * x50 + x11 * x36 * x9 + x15 * x37 + x21 * x49 + x29 * x5 - x33 * x5;
                T x52 = x42 * x51;
                T x53 = x0 * x24;
                T x54 = x19 * x37;
                T x55 = x0 * x21;
                T x56 = x1 * x53 + x11 * x54 - x19 * x55 + x31 * x5 - x34 * x5 + x36 * x39;
                T x57 = x42 * x56;
                T x58 = x18 * x30;
                T x59 = x19 * x50 - x24 * x49 + x28 * x5 + x30 * x37 * x9 - x5 * x58 - x54 * x7;
                T x60 = x42 * x59;
                T x61 = x40 * x42;
                T x62 = -x52 - x57 - x60 - x61;
                T x63 = -x10 * x28 + x10 * x58 + x2 * x30 * x36 + x2 * x53 + x22 * x37 * x7 - x22 * x50;
                T x64 = x42 * x63;
                T x65 = x35 * x42;
                T x66 = -x10 * x29 + x10 * x33 + x13 * x37 - x2 * x55 - x37 * x8 + x50 * x6;
                T x67 = x42 * x66;
                T x68 = -x10 * x31 + x10 * x34 + x22 * x55 + x32 * x37 - x37 * x44 - x53 * x6;
                T x69 = x42 * x68;
                T x70 = -x64 - x65 - x67 - x69;
                T x71 = x18 * x5;
                T x72 = x10 * x9;
                T x73 = x10 * x18;
                T x74 = x2 * x5;
                T x75 = -x19 * x73 + x22 * x71 + x23 * x36 + x24 * x72 - x24 * x74 + x26 * x37;
                T x76 = x42 * x75;
                T x77 = x27 * x42;
                T x78 = x1 * x73 + x14 * x37 - x21 * x72 + x21 * x74 - x3 * x37 - x6 * x71;
                T x79 = x42 * x78;
                T x80 = -x1 * x10 * x24 + x10 * x19 * x21 - x20 * x37 - x21 * x22 * x5 + x24 * x5 * x6 + x25 * x37;
                T x81 = x42 * x80;
                T x82 = -x76 - x77 - x79 - x81;
                T x83 = weight * x41;
                T x84 = x83 * (x47 * x48 + x61 * x62 + x65 * x70 + x77 * x82);
                T x85 = x83 * (x46 * x48 + x57 * x62 + x69 * x70 + x81 * x82);
                T x86 = x83 * (x60 * x62 + x64 * x70 + x76 * x82);
                T x87 = x83 * (x43 * x48 + x52 * x62 + x67 * x70 + x79 * x82);
                T x88 = pow(x41, -2);
                T x89 = x40 * x88;
                T x90 = x35 * x88;
                T x91 = x38 * x88;
                T x92 = x27 * x88;
                T x93 = x83 * (x45 * x91 + x56 * x89 + x68 * x90 + x80 * x92);
                T x94 = x83 * (x59 * x89 + x63 * x90 + x75 * x92);
                T x95 = x83 * (x16 * x91 + x51 * x89 + x66 * x90 + x78 * x92);
                T x96 = x56 * x88;
                T x97 = x68 * x88;
                T x98 = x80 * x88;
                T x99 = x83 * (x59 * x96 + x63 * x97 + x75 * x98);
                T x100 = x83 * (x16 * x45 * x88 + x51 * x96 + x66 * x97 + x78 * x98);
                T x101 = x83 * (x51 * x59 * x88 + x63 * x66 * x88 + x75 * x78 * x88);
                T x102 = u[0] * x48 + u[1] * x47 + u[2] * x46 + u[4] * x43;
                T x103 = u[0] * x62 + u[1] * x61 + u[2] * x57 + u[3] * x60 + u[4] * x52;
                T x104 = u[0] * x70 + u[1] * x65 + u[2] * x69 + u[3] * x64 + u[4] * x67;
                T x105 = u[0] * x82 + u[1] * x77 + u[2] * x81 + u[3] * x76 + u[4] * x79;
                T x106 = 5 * x83;
                H[0] += x83 * (pow(x48, 2) + pow(x62, 2) + pow(x70, 2) + pow(x82, 2));
                H[1] += x84;
                H[2] += x85;
                H[3] += x86;
                H[4] += x87;
                H[5] += x84;
                H[6] += x83 * (pow(x27, 2) * x88 + pow(x35, 2) * x88 + pow(x38, 2) * x88 + pow(x40, 2) * x88);
                H[7] += x93;
                H[8] += x94;
                H[9] += x95;
                H[10] += x85;
                H[11] += x93;
                H[12] += x83 * (pow(x45, 2) * x88 + pow(x56, 2) * x88 + pow(x68, 2) * x88 + pow(x80, 2) * x88);
                H[13] += x99;
                H[14] += x100;
                H[15] += x86;
                H[16] += x94;
                H[17] += x99;
                H[18] += x83 * (pow(x59, 2) * x88 + pow(x63, 2) * x88 + pow(x75, 2) * x88);
                H[19] += x101;
                H[20] += x87;
                H[21] += x95;
                H[22] += x100;
                H[23] += x101;
                H[24] += x83 * (pow(x16, 2) * x88 + pow(x51, 2) * x88 + pow(x66, 2) * x88 + pow(x78, 2) * x88);
                g[0] += x106 * (x102 * x48 + x103 * x62 + x104 * x70 + x105 * x82);
                g[1] += x106 * (x102 * x47 + x103 * x61 + x104 * x65 + x105 * x77);
                g[2] += x106 * (x102 * x46 + x103 * x57 + x104 * x69 + x105 * x81);
                g[3] += x106 * (x103 * x60 + x104 * x64 + x105 * x76);
                g[4] += x106 * (x102 * x43 + x103 * x52 + x104 * x67 + x105 * x79);
                e += x83 * (pow(x102, 2) + pow(x103, 2) + pow(x104, 2) + pow(x105, 2));
            }

            // TODO
        };
    }  // namespace kernels

    namespace kokkos {
        template <class FunctionSpace, class FE>
        using LaplaceOperatorPentatope5 = utopia::kokkos::AutoKernel<FunctionSpace,
            FE,
            utopia::kernels::LaplaceOperator<utopia::kernels::Pentatope5<typename FE::Scalar, typename FE::Scalar>>,
            4>;
    }
}  // namespace utopia

#endif  // UTOPIA_TPL_MATERIAL_LaplaceOperator_Pentatope5_4_IMPL_hpp
