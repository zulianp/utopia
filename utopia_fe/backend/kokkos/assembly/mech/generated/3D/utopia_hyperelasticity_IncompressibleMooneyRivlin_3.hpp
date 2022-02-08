#ifndef UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_3_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_3_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_IncompressibleMooneyRivlin.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of IncompressibleMooneyRivlin for dimension 3
         */
        template <typename T>
        class IncompressibleMooneyRivlin<T, 3> {
        public:
            static constexpr int Dim = 3;

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    in.get("C1", C1);
                    in.get("C2", C2);
                }

                T C1{1.0};
                T C2{1.0};
            };

            IncompressibleMooneyRivlin(const Params &params) {
                C1 = params.C1;
                C2 = params.C2;
            }

            UTOPIA_FUNCTION void hessian(const T *UTOPIA_RESTRICT f,
                                         const T p,
                                         const T *grad_test,
                                         const T *grad_trial,
                                         const T fun_test,
                                         const T fun_trial,
                                         const T dx,
                                         T *UTOPIA_RESTRICT bf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = 2 * f[3];
                T x1 = -f[4] * x0;
                T x2 = 2 * f[6];
                T x3 = -f[7] * x2;
                T x4 = C2 * (x1 + x3);
                T x5 = -f[5] * x0;
                T x6 = -f[8] * x2;
                T x7 = x5 + x6;
                T x8 = C2 * grad_test[2];
                T x9 = 2 * C1;
                T x10 = 2 * pow(f[5], 2);
                T x11 = 2 * pow(f[8], 2);
                T x12 = x10 + x11;
                T x13 = 2 * pow(f[4], 2);
                T x14 = 2 * pow(f[7], 2);
                T x15 = x13 + x14;
                T x16 = 2 * f[4];
                T x17 = -f[5] * x16;
                T x18 = 2 * f[7];
                T x19 = -f[8] * x18;
                T x20 = x17 + x19;
                T x21 = 2 * pow(f[3], 2);
                T x22 = 2 * pow(f[6], 2);
                T x23 = x21 + x22;
                T x24 = C2 * grad_test[0];
                T x25 = C2 * grad_test[1];
                T x26 = -f[1] * x16;
                T x27 = 2 * f[2];
                T x28 = -f[5] * x27;
                T x29 = x24 * (x26 + x28);
                T x30 = f[7] * p;
                T x31 = f[0] * f[5];
                T x32 = f[2] * f[3];
                T x33 = C2 * (-2 * x31 + 4 * x32) + x30;
                T x34 = f[8] * p;
                T x35 = f[0] * f[4];
                T x36 = f[1] * f[3];
                T x37 = C2 * (-2 * x35 + 4 * x36) - x34;
                T x38 = -f[0] * x0;
                T x39 = x25 * (x28 + x38);
                T x40 = C2 * (4 * x35 - 2 * x36) + x34;
                T x41 = f[6] * p;
                T x42 = f[1] * f[5];
                T x43 = f[2] * f[4];
                T x44 = C2 * (-2 * x42 + 4 * x43) - x41;
                T x45 = x8 * (x26 + x38);
                T x46 = C2 * (4 * x42 - 2 * x43) + x41;
                T x47 = C2 * (4 * x31 - 2 * x32) - x30;
                T x48 = -f[1] * x18;
                T x49 = -f[8] * x27;
                T x50 = x24 * (x48 + x49);
                T x51 = f[5] * p;
                T x52 = f[0] * f[7];
                T x53 = f[1] * f[6];
                T x54 = C2 * (-2 * x52 + 4 * x53) + x51;
                T x55 = f[4] * p;
                T x56 = f[0] * f[8];
                T x57 = f[2] * f[6];
                T x58 = C2 * (-2 * x56 + 4 * x57) - x55;
                T x59 = -f[0] * x2;
                T x60 = x25 * (x49 + x59);
                T x61 = f[3] * p;
                T x62 = f[1] * f[8];
                T x63 = f[2] * f[7];
                T x64 = C2 * (-2 * x62 + 4 * x63) + x61;
                T x65 = C2 * (4 * x52 - 2 * x53) - x51;
                T x66 = x8 * (x48 + x59);
                T x67 = C2 * (4 * x56 - 2 * x57) + x55;
                T x68 = C2 * (4 * x62 - 2 * x63) - x61;
                T x69 = -2 * f[0] * f[1];
                T x70 = x3 + x69;
                T x71 = -f[0] * x27;
                T x72 = x6 + x71;
                T x73 = 2 * pow(f[2], 2);
                T x74 = x11 + x73;
                T x75 = 2 * pow(f[1], 2);
                T x76 = x14 + x75;
                T x77 = -f[1] * x27;
                T x78 = x19 + x77;
                T x79 = 2 * pow(f[0], 2);
                T x80 = x22 + x79;
                T x81 = -f[7] * x16;
                T x82 = -2 * f[5] * f[8];
                T x83 = x24 * (x81 + x82);
                T x84 = f[1] * p;
                T x85 = f[3] * f[8];
                T x86 = f[5] * f[6];
                T x87 = C2 * (-2 * x85 + 4 * x86) + x84;
                T x88 = f[2] * p;
                T x89 = f[3] * f[7];
                T x90 = f[4] * f[6];
                T x91 = C2 * (-2 * x89 + 4 * x90) - x88;
                T x92 = -f[6] * x0;
                T x93 = x25 * (x82 + x92);
                T x94 = C2 * (4 * x89 - 2 * x90) + x88;
                T x95 = f[0] * p;
                T x96 = f[4] * f[8];
                T x97 = f[5] * f[7];
                T x98 = C2 * (-2 * x96 + 4 * x97) - x95;
                T x99 = x8 * (x81 + x92);
                T x100 = C2 * (4 * x96 - 2 * x97) + x95;
                T x101 = C2 * (4 * x85 - 2 * x86) - x84;
                T x102 = x1 + x69;
                T x103 = x5 + x71;
                T x104 = x10 + x73;
                T x105 = x13 + x75;
                T x106 = x17 + x77;
                T x107 = x21 + x79;
                T x108 = x96 - x97;
                T x109 = x89 - x90;
                T x110 = x85 - x86;
                T x111 = dx * fun_test;
                T x112 = x56 - x57;
                T x113 = x62 - x63;
                T x114 = x52 - x53;
                T x115 = x42 - x43;
                T x116 = x35 - x36;
                T x117 = x31 - x32;
                T x118 = dx * fun_trial;
                bf[0] += dx * (grad_trial[0] * (grad_test[0] * (C2 * (x12 + x15) + x9) + grad_test[1] * x4 + x7 * x8) +
                               grad_trial[1] * (grad_test[0] * x4 + grad_test[1] * (C2 * (x12 + x23) + x9) + x20 * x8) +
                               grad_trial[2] * (grad_test[2] * (C2 * (x15 + x23) + x9) + x20 * x25 + x24 * x7));
                bf[1] += dx * (grad_trial[0] * (grad_test[1] * x37 + grad_test[2] * x33 + x29) +
                               grad_trial[1] * (grad_test[0] * x40 + grad_test[2] * x44 + x39) +
                               grad_trial[2] * (grad_test[0] * x47 + grad_test[1] * x46 + x45));
                bf[2] += dx * (grad_trial[0] * (grad_test[1] * x54 + grad_test[2] * x58 + x50) +
                               grad_trial[1] * (grad_test[0] * x65 + grad_test[2] * x64 + x60) +
                               grad_trial[2] * (grad_test[0] * x67 + grad_test[1] * x68 + x66));
                bf[4] += dx * (grad_trial[0] * (grad_test[1] * x40 + grad_test[2] * x47 + x29) +
                               grad_trial[1] * (grad_test[0] * x37 + grad_test[2] * x46 + x39) +
                               grad_trial[2] * (grad_test[0] * x33 + grad_test[1] * x44 + x45));
                bf[5] += dx * (grad_trial[0] * (grad_test[0] * (C2 * (x74 + x76) + x9) + x25 * x70 + x72 * x8) +
                               grad_trial[1] * (grad_test[1] * (C2 * (x74 + x80) + x9) + x24 * x70 + x78 * x8) +
                               grad_trial[2] * (grad_test[2] * (C2 * (x76 + x80) + x9) + x24 * x72 + x25 * x78));
                bf[6] += dx * (grad_trial[0] * (grad_test[1] * x91 + grad_test[2] * x87 + x83) +
                               grad_trial[1] * (grad_test[0] * x94 + grad_test[2] * x98 + x93) +
                               grad_trial[2] * (grad_test[0] * x101 + grad_test[1] * x100 + x99));
                bf[8] += dx * (grad_trial[0] * (grad_test[1] * x65 + grad_test[2] * x67 + x50) +
                               grad_trial[1] * (grad_test[0] * x54 + grad_test[2] * x68 + x60) +
                               grad_trial[2] * (grad_test[0] * x58 + grad_test[1] * x64 + x66));
                bf[9] += dx * (grad_trial[0] * (grad_test[1] * x94 + grad_test[2] * x101 + x83) +
                               grad_trial[1] * (grad_test[0] * x91 + grad_test[2] * x100 + x93) +
                               grad_trial[2] * (grad_test[0] * x87 + grad_test[1] * x98 + x99));
                bf[10] += dx * (grad_trial[0] * (grad_test[0] * (C2 * (x104 + x105) + x9) + x102 * x25 + x103 * x8) +
                                grad_trial[1] * (grad_test[1] * (C2 * (x104 + x107) + x9) + x102 * x24 + x106 * x8) +
                                grad_trial[2] * (grad_test[2] * (C2 * (x105 + x107) + x9) + x103 * x24 + x106 * x25));
                bf[15] += 0;
                bf[12] += x111 * (grad_trial[0] * x108 - grad_trial[1] * x110 + grad_trial[2] * x109);
                bf[13] += x111 * (-grad_trial[0] * x113 + grad_trial[1] * x112 - grad_trial[2] * x114);
                bf[14] += x111 * (grad_trial[0] * x115 - grad_trial[1] * x117 + grad_trial[2] * x116);
                bf[3] += x118 * (grad_test[0] * x108 - grad_test[1] * x110 + grad_test[2] * x109);
                bf[7] += x118 * (-grad_test[0] * x113 + grad_test[1] * x112 - grad_test[2] * x114);
                bf[11] += x118 * (grad_test[0] * x115 - grad_test[1] * x117 + grad_test[2] * x116);
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T p,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T fun_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = 2 * f[0];
                T x1 = f[4] * f[8];
                T x2 = f[5] * f[7];
                T x3 = pow(f[0], 2) + pow(f[3], 2) + pow(f[6], 2);
                T x4 = f[0] * f[1] + f[3] * f[4] + f[6] * f[7];
                T x5 = 2 * f[1];
                T x6 = f[0] * f[2] + f[3] * f[5] + f[6] * f[8];
                T x7 = 2 * f[2];
                T x8 = pow(f[1], 2) + pow(f[4], 2) + pow(f[7], 2);
                T x9 = pow(f[2], 2) + pow(f[5], 2) + pow(f[8], 2);
                T x10 = x3 + x8 + x9;
                T x11 = f[5] * f[6];
                T x12 = f[3] * f[8];
                T x13 = f[1] * f[2] + f[4] * f[5] + f[7] * f[8];
                T x14 = f[3] * f[7];
                T x15 = f[4] * f[6];
                T x16 = 2 * f[3];
                T x17 = 2 * f[4];
                T x18 = 2 * f[5];
                T x19 = 2 * f[6];
                T x20 = 2 * f[7];
                T x21 = 2 * f[8];
                lf[0] +=
                    dx * (grad_test[0] * (C1 * x0 + C2 * (x0 * x10 - x0 * x3 - x4 * x5 - x6 * x7) + p * (x1 - x2)) +
                          grad_test[1] * (C1 * x5 + C2 * (-x0 * x4 + x10 * x5 - x13 * x7 - x5 * x8) + p * (x11 - x12)) +
                          grad_test[2] * (C1 * x7 + C2 * (-x0 * x6 + x10 * x7 - x13 * x5 - x7 * x9) + p * (x14 - x15)));
                lf[1] += dx * (grad_test[0] * (C1 * x16 + C2 * (x10 * x16 - x16 * x3 - x17 * x4 - x18 * x6) +
                                               p * (-f[1] * f[8] + f[2] * f[7])) +
                               grad_test[1] * (C1 * x17 + C2 * (x10 * x17 - x13 * x18 - x16 * x4 - x17 * x8) +
                                               p * (f[0] * f[8] - f[2] * f[6])) +
                               grad_test[2] * (C1 * x18 + C2 * (x10 * x18 - x13 * x17 - x16 * x6 - x18 * x9) +
                                               p * (-f[0] * f[7] + f[1] * f[6])));
                lf[2] += dx * (grad_test[0] * (C1 * x19 + C2 * (x10 * x19 - x19 * x3 - x20 * x4 - x21 * x6) +
                                               p * (f[1] * f[5] - f[2] * f[4])) +
                               grad_test[1] * (C1 * x20 + C2 * (x10 * x20 - x13 * x21 - x19 * x4 - x20 * x8) +
                                               p * (-f[0] * f[5] + f[2] * f[3])) +
                               grad_test[2] * (C1 * x21 + C2 * (x10 * x21 - x13 * x20 - x19 * x6 - x21 * x9) +
                                               p * (f[0] * f[4] - f[1] * f[3])));
                lf[3] +=
                    dx * fun_test * (f[0] * x1 - f[0] * x2 + f[1] * x11 - f[1] * x12 + f[2] * x14 - f[2] * x15 - 1);
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T p, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = pow(f[0], 2) + pow(f[3], 2) + pow(f[6], 2);
                T x1 = pow(f[1], 2) + pow(f[4], 2) + pow(f[7], 2);
                T x2 = pow(f[2], 2) + pow(f[5], 2) + pow(f[8], 2);
                T x3 = x0 + x1 + x2;
                e += dx * (C1 * (x3 - 3) +
                           C2 * (-1.0 / 2.0 * pow(x0, 2) - 1.0 / 2.0 * pow(x1, 2) - 1.0 / 2.0 * pow(x2, 2) +
                                 (1.0 / 2.0) * pow(x3, 2) - pow(f[0] * f[1] + f[3] * f[4] + f[6] * f[7], 2) -
                                 pow(f[0] * f[2] + f[3] * f[5] + f[6] * f[8], 2) -
                                 pow(f[1] * f[2] + f[4] * f[5] + f[7] * f[8], 2) - 3) +
                           p * (f[0] * f[4] * f[8] - f[0] * f[5] * f[7] - f[1] * f[3] * f[8] + f[1] * f[5] * f[6] +
                                f[2] * f[3] * f[7] - f[2] * f[4] * f[6] - 1));
            }

            UTOPIA_FUNCTION void eval(const T *UTOPIA_RESTRICT f,
                                      const T p,
                                      const T *grad_test,
                                      const T *grad_trial,
                                      const T fun_test,
                                      const T fun_trial,
                                      const T dx,
                                      T &e,
                                      T *UTOPIA_RESTRICT lf,
                                      T *UTOPIA_RESTRICT bf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = pow(f[0], 2);
                T x1 = pow(f[3], 2);
                T x2 = pow(f[6], 2);
                T x3 = x0 + x1 + x2;
                T x4 = pow(f[1], 2);
                T x5 = pow(f[4], 2);
                T x6 = pow(f[7], 2);
                T x7 = x4 + x5 + x6;
                T x8 = pow(f[2], 2);
                T x9 = pow(f[5], 2);
                T x10 = pow(f[8], 2);
                T x11 = x10 + x8 + x9;
                T x12 = x11 + x3 + x7;
                T x13 = f[4] * f[8];
                T x14 = f[5] * f[6];
                T x15 = f[3] * f[7];
                T x16 = f[5] * f[7];
                T x17 = f[3] * f[8];
                T x18 = f[4] * f[6];
                T x19 = f[0] * x13 - f[0] * x16 + f[1] * x14 - f[1] * x17 + f[2] * x15 - f[2] * x18 - 1;
                T x20 = f[0] * f[1];
                T x21 = f[3] * f[4];
                T x22 = f[6] * f[7];
                T x23 = x20 + x21 + x22;
                T x24 = f[0] * f[2];
                T x25 = f[3] * f[5];
                T x26 = f[6] * f[8];
                T x27 = x24 + x25 + x26;
                T x28 = f[1] * f[2];
                T x29 = f[4] * f[5];
                T x30 = f[7] * f[8];
                T x31 = x28 + x29 + x30;
                T x32 = 2 * C1;
                T x33 = x13 - x16;
                T x34 = 2 * f[0];
                T x35 = 2 * f[1];
                T x36 = 2 * f[2];
                T x37 = x15 - x18;
                T x38 = f[2] * f[7];
                T x39 = f[1] * f[8];
                T x40 = 2 * f[3];
                T x41 = 2 * f[4];
                T x42 = 2 * f[5];
                T x43 = f[0] * f[8];
                T x44 = f[2] * f[6];
                T x45 = x43 - x44;
                T x46 = f[1] * f[6];
                T x47 = f[0] * f[7];
                T x48 = f[1] * f[5];
                T x49 = f[2] * f[4];
                T x50 = x48 - x49;
                T x51 = 2 * f[6];
                T x52 = 2 * f[7];
                T x53 = 2 * f[8];
                T x54 = f[2] * f[3];
                T x55 = f[0] * f[5];
                T x56 = f[0] * f[4];
                T x57 = f[1] * f[3];
                T x58 = x56 - x57;
                T x59 = dx * fun_test;
                T x60 = -2 * x21;
                T x61 = -2 * x22;
                T x62 = C2 * (x60 + x61);
                T x63 = -2 * x25;
                T x64 = -2 * x26;
                T x65 = x63 + x64;
                T x66 = C2 * grad_test[2];
                T x67 = 2 * x9;
                T x68 = 2 * x10;
                T x69 = x67 + x68;
                T x70 = 2 * x5;
                T x71 = 2 * x6;
                T x72 = x70 + x71;
                T x73 = -2 * x29;
                T x74 = -2 * x30;
                T x75 = x73 + x74;
                T x76 = 2 * x1;
                T x77 = 2 * x2;
                T x78 = x76 + x77;
                T x79 = C2 * grad_test[0];
                T x80 = C2 * grad_test[1];
                T x81 = -f[4] * x35;
                T x82 = -f[5] * x36;
                T x83 = x79 * (x81 + x82);
                T x84 = f[7] * p;
                T x85 = C2 * (4 * x54 - 2 * x55) + x84;
                T x86 = f[8] * p;
                T x87 = C2 * (-2 * x56 + 4 * x57) - x86;
                T x88 = -f[3] * x34;
                T x89 = x80 * (x82 + x88);
                T x90 = C2 * (4 * x56 - 2 * x57) + x86;
                T x91 = f[6] * p;
                T x92 = C2 * (-2 * x48 + 4 * x49) - x91;
                T x93 = x66 * (x81 + x88);
                T x94 = C2 * (4 * x48 - 2 * x49) + x91;
                T x95 = C2 * (-2 * x54 + 4 * x55) - x84;
                T x96 = -f[7] * x35;
                T x97 = -f[8] * x36;
                T x98 = x79 * (x96 + x97);
                T x99 = f[5] * p;
                T x100 = C2 * (4 * x46 - 2 * x47) + x99;
                T x101 = f[4] * p;
                T x102 = C2 * (-2 * x43 + 4 * x44) - x101;
                T x103 = -f[6] * x34;
                T x104 = x80 * (x103 + x97);
                T x105 = f[3] * p;
                T x106 = C2 * (4 * x38 - 2 * x39) + x105;
                T x107 = C2 * (-2 * x46 + 4 * x47) - x99;
                T x108 = x66 * (x103 + x96);
                T x109 = C2 * (4 * x43 - 2 * x44) + x101;
                T x110 = C2 * (-2 * x38 + 4 * x39) - x105;
                T x111 = -2 * x20;
                T x112 = x111 + x61;
                T x113 = -2 * x24;
                T x114 = x113 + x64;
                T x115 = 2 * x8;
                T x116 = x115 + x68;
                T x117 = 2 * x4;
                T x118 = x117 + x71;
                T x119 = -2 * x28;
                T x120 = x119 + x74;
                T x121 = 2 * x0;
                T x122 = x121 + x77;
                T x123 = -f[7] * x41;
                T x124 = -f[8] * x42;
                T x125 = x79 * (x123 + x124);
                T x126 = f[1] * p;
                T x127 = C2 * (4 * x14 - 2 * x17) + x126;
                T x128 = f[2] * p;
                T x129 = C2 * (-2 * x15 + 4 * x18) - x128;
                T x130 = -f[6] * x40;
                T x131 = x80 * (x124 + x130);
                T x132 = C2 * (4 * x15 - 2 * x18) + x128;
                T x133 = f[0] * p;
                T x134 = C2 * (-2 * x13 + 4 * x16) - x133;
                T x135 = x66 * (x123 + x130);
                T x136 = C2 * (4 * x13 - 2 * x16) + x133;
                T x137 = C2 * (-2 * x14 + 4 * x17) - x126;
                T x138 = x111 + x60;
                T x139 = x113 + x63;
                T x140 = x115 + x67;
                T x141 = x117 + x70;
                T x142 = x119 + x73;
                T x143 = x121 + x76;
                T x144 = -x14 + x17;
                T x145 = -x38 + x39;
                T x146 = -x46 + x47;
                T x147 = -x54 + x55;
                T x148 = dx * fun_trial;
                e += dx * (C1 * (x12 - 3) +
                           C2 * (-1.0 / 2.0 * pow(x11, 2) + (1.0 / 2.0) * pow(x12, 2) - pow(x23, 2) - pow(x27, 2) -
                                 1.0 / 2.0 * pow(x3, 2) - pow(x31, 2) - 1.0 / 2.0 * pow(x7, 2) - 3) +
                           p * x19);
                lf[0] +=
                    dx *
                    (grad_test[0] * (C2 * (x12 * x34 - x23 * x35 - x27 * x36 - x3 * x34) + f[0] * x32 + p * x33) +
                     grad_test[1] *
                         (C2 * (x12 * x35 - x23 * x34 - x31 * x36 - x35 * x7) + f[1] * x32 + p * (x14 - x17)) +
                     grad_test[2] * (C2 * (-x11 * x36 + x12 * x36 - x27 * x34 - x31 * x35) + f[2] * x32 + p * x37));
                lf[1] +=
                    dx * (grad_test[0] *
                              (C2 * (x12 * x40 - x23 * x41 - x27 * x42 - x3 * x40) + f[3] * x32 + p * (x38 - x39)) +
                          grad_test[1] * (C2 * (x12 * x41 - x23 * x40 - x31 * x42 - x41 * x7) + f[4] * x32 + p * x45) +
                          grad_test[2] *
                              (C2 * (-x11 * x42 + x12 * x42 - x27 * x40 - x31 * x41) + f[5] * x32 + p * (x46 - x47)));
                lf[2] +=
                    dx *
                    (grad_test[0] * (C2 * (x12 * x51 - x23 * x52 - x27 * x53 - x3 * x51) + f[6] * x32 + p * x50) +
                     grad_test[1] *
                         (C2 * (x12 * x52 - x23 * x51 - x31 * x53 - x52 * x7) + f[7] * x32 + p * (x54 - x55)) +
                     grad_test[2] * (C2 * (-x11 * x53 + x12 * x53 - x27 * x51 - x31 * x52) + f[8] * x32 + p * x58));
                lf[3] += x19 * x59;
                bf[0] +=
                    dx * (grad_trial[0] * (grad_test[0] * (C2 * (x69 + x72) + x32) + grad_test[1] * x62 + x65 * x66) +
                          grad_trial[1] * (grad_test[0] * x62 + grad_test[1] * (C2 * (x69 + x78) + x32) + x66 * x75) +
                          grad_trial[2] * (grad_test[2] * (C2 * (x72 + x78) + x32) + x65 * x79 + x75 * x80));
                bf[1] += dx * (grad_trial[0] * (grad_test[1] * x87 + grad_test[2] * x85 + x83) +
                               grad_trial[1] * (grad_test[0] * x90 + grad_test[2] * x92 + x89) +
                               grad_trial[2] * (grad_test[0] * x95 + grad_test[1] * x94 + x93));
                bf[2] += dx * (grad_trial[0] * (grad_test[1] * x100 + grad_test[2] * x102 + x98) +
                               grad_trial[1] * (grad_test[0] * x107 + grad_test[2] * x106 + x104) +
                               grad_trial[2] * (grad_test[0] * x109 + grad_test[1] * x110 + x108));
                bf[4] += dx * (grad_trial[0] * (grad_test[1] * x90 + grad_test[2] * x95 + x83) +
                               grad_trial[1] * (grad_test[0] * x87 + grad_test[2] * x94 + x89) +
                               grad_trial[2] * (grad_test[0] * x85 + grad_test[1] * x92 + x93));
                bf[5] += dx * (grad_trial[0] * (grad_test[0] * (C2 * (x116 + x118) + x32) + x112 * x80 + x114 * x66) +
                               grad_trial[1] * (grad_test[1] * (C2 * (x116 + x122) + x32) + x112 * x79 + x120 * x66) +
                               grad_trial[2] * (grad_test[2] * (C2 * (x118 + x122) + x32) + x114 * x79 + x120 * x80));
                bf[6] += dx * (grad_trial[0] * (grad_test[1] * x129 + grad_test[2] * x127 + x125) +
                               grad_trial[1] * (grad_test[0] * x132 + grad_test[2] * x134 + x131) +
                               grad_trial[2] * (grad_test[0] * x137 + grad_test[1] * x136 + x135));
                bf[8] += dx * (grad_trial[0] * (grad_test[1] * x107 + grad_test[2] * x109 + x98) +
                               grad_trial[1] * (grad_test[0] * x100 + grad_test[2] * x110 + x104) +
                               grad_trial[2] * (grad_test[0] * x102 + grad_test[1] * x106 + x108));
                bf[9] += dx * (grad_trial[0] * (grad_test[1] * x132 + grad_test[2] * x137 + x125) +
                               grad_trial[1] * (grad_test[0] * x129 + grad_test[2] * x136 + x131) +
                               grad_trial[2] * (grad_test[0] * x127 + grad_test[1] * x134 + x135));
                bf[10] += dx * (grad_trial[0] * (grad_test[0] * (C2 * (x140 + x141) + x32) + x138 * x80 + x139 * x66) +
                                grad_trial[1] * (grad_test[1] * (C2 * (x140 + x143) + x32) + x138 * x79 + x142 * x66) +
                                grad_trial[2] * (grad_test[2] * (C2 * (x141 + x143) + x32) + x139 * x79 + x142 * x80));
                bf[15] += 0;
                bf[12] += x59 * (grad_trial[0] * x33 - grad_trial[1] * x144 + grad_trial[2] * x37);
                bf[13] += x59 * (-grad_trial[0] * x145 + grad_trial[1] * x45 - grad_trial[2] * x146);
                bf[14] += x59 * (grad_trial[0] * x50 - grad_trial[1] * x147 + grad_trial[2] * x58);
                bf[3] += x148 * (grad_test[0] * x33 - grad_test[1] * x144 + grad_test[2] * x37);
                bf[7] += x148 * (-grad_test[0] * x145 + grad_test[1] * x45 - grad_test[2] * x146);
                bf[11] += x148 * (grad_test[0] * x50 - grad_test[1] * x147 + grad_test[2] * x58);
            }

            T C1{1.0};
            T C2{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_3_IMPL_hpp
