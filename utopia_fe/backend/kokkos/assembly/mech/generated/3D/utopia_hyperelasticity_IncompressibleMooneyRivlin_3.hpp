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
                T x0 = f[3] * f[4];
                T x1 = f[6] * f[7];
                T x2 = C2 * (x0 + x1);
                T x3 = f[3] * f[5];
                T x4 = f[6] * f[8];
                T x5 = x3 + x4;
                T x6 = C2 * grad_test[2];
                T x7 = pow(f[5], 2);
                T x8 = pow(f[8], 2);
                T x9 = x7 + x8;
                T x10 = pow(f[4], 2);
                T x11 = pow(f[7], 2);
                T x12 = x10 + x11;
                T x13 = f[4] * f[5];
                T x14 = f[7] * f[8];
                T x15 = x13 + x14;
                T x16 = pow(f[3], 2);
                T x17 = pow(f[6], 2);
                T x18 = x16 + x17;
                T x19 = C2 * grad_test[0];
                T x20 = C2 * grad_test[1];
                T x21 = 2 * dx;
                T x22 = f[1] * f[4];
                T x23 = f[2] * f[5];
                T x24 = 2 * x19;
                T x25 = x24 * (x22 + x23);
                T x26 = f[8] * p;
                T x27 = f[0] * f[4];
                T x28 = f[1] * f[3];
                T x29 = 2 * C2;
                T x30 = x26 + x29 * (x27 - 2 * x28);
                T x31 = f[7] * p;
                T x32 = f[0] * f[5];
                T x33 = f[2] * f[3];
                T x34 = -x29 * (x32 - 2 * x33) + x31;
                T x35 = f[0] * f[3];
                T x36 = 2 * x20;
                T x37 = x36 * (x23 + x35);
                T x38 = f[6] * p;
                T x39 = f[1] * f[5];
                T x40 = f[2] * f[4];
                T x41 = x29 * (x39 - 2 * x40) + x38;
                T x42 = -x28;
                T x43 = x26 + x29 * (2 * x27 + x42);
                T x44 = 2 * x6;
                T x45 = x44 * (x22 + x35);
                T x46 = -x40;
                T x47 = x29 * (2 * x39 + x46) + x38;
                T x48 = -x33;
                T x49 = x29 * (2 * x32 + x48) - x31;
                T x50 = f[1] * f[7];
                T x51 = f[2] * f[8];
                T x52 = x24 * (x50 + x51);
                T x53 = f[4] * p;
                T x54 = f[0] * f[8];
                T x55 = f[2] * f[6];
                T x56 = x29 * (x54 - 2 * x55) + x53;
                T x57 = f[5] * p;
                T x58 = f[0] * f[7];
                T x59 = f[1] * f[6];
                T x60 = x29 * (x58 - 2 * x59);
                T x61 = f[0] * f[6];
                T x62 = x44 * (x50 + x61);
                T x63 = -x55;
                T x64 = x29 * (2 * x54 + x63) + x53;
                T x65 = f[3] * p;
                T x66 = -x65;
                T x67 = f[1] * f[8];
                T x68 = f[2] * f[7];
                T x69 = -x68;
                T x70 = x29 * (2 * x67 + x69) + x66;
                T x71 = -x36 * (x51 + x61);
                T x72 = -x57;
                T x73 = -x59;
                T x74 = x29 * (2 * x58 + x73) + x72;
                T x75 = x29 * (x67 - 2 * x68);
                T x76 = f[0] * f[1];
                T x77 = x1 + x76;
                T x78 = f[0] * f[2];
                T x79 = x4 + x78;
                T x80 = pow(f[2], 2);
                T x81 = x8 + x80;
                T x82 = pow(f[1], 2);
                T x83 = x11 + x82;
                T x84 = f[1] * f[2];
                T x85 = x14 + x84;
                T x86 = pow(f[0], 2);
                T x87 = x17 + x86;
                T x88 = f[4] * f[7];
                T x89 = f[5] * f[8];
                T x90 = x24 * (x88 + x89);
                T x91 = f[2] * p;
                T x92 = f[3] * f[7];
                T x93 = f[4] * f[6];
                T x94 = x29 * (x92 - 2 * x93) + x91;
                T x95 = f[1] * p;
                T x96 = f[3] * f[8];
                T x97 = f[5] * f[6];
                T x98 = -x29 * (x96 - 2 * x97) + x95;
                T x99 = f[3] * f[6];
                T x100 = x36 * (x89 + x99);
                T x101 = f[0] * p;
                T x102 = f[4] * f[8];
                T x103 = f[5] * f[7];
                T x104 = x101 + x29 * (x102 - 2 * x103);
                T x105 = -x93;
                T x106 = x29 * (x105 + 2 * x92) + x91;
                T x107 = x44 * (x88 + x99);
                T x108 = -x103;
                T x109 = x101 + x29 * (2 * x102 + x108);
                T x110 = -x97;
                T x111 = x29 * (x110 + 2 * x96) - x95;
                T x112 = x0 + x76;
                T x113 = x3 + x78;
                T x114 = x7 + x80;
                T x115 = x10 + x82;
                T x116 = x13 + x84;
                T x117 = x16 + x86;
                T x118 = x102 + x108;
                T x119 = x105 + x92;
                T x120 = x110 + x96;
                T x121 = dx * fun_test;
                T x122 = x67 + x69;
                T x123 = x58 + x73;
                T x124 = x54 + x63;
                T x125 = x39 + x46;
                T x126 = x27 + x42;
                T x127 = x32 + x48;
                T x128 = dx * fun_trial;
                bf[0] +=
                    -x21 * (grad_trial[0] * (-grad_test[0] * (C1 + C2 * (x12 + x9)) + grad_test[1] * x2 + x5 * x6) +
                            grad_trial[1] * (grad_test[0] * x2 - grad_test[1] * (C1 + C2 * (x18 + x9)) + x15 * x6) +
                            grad_trial[2] * (-grad_test[2] * (C1 + C2 * (x12 + x18)) + x15 * x20 + x19 * x5));
                bf[1] += dx * (-grad_trial[0] * (grad_test[1] * x30 - grad_test[2] * x34 + x25) -
                               grad_trial[1] * (-grad_test[0] * x43 + grad_test[2] * x41 + x37) +
                               grad_trial[2] * (grad_test[0] * x49 + grad_test[1] * x47 - x45));
                bf[2] += dx * (-grad_trial[0] * (-grad_test[1] * (x57 - x60) + grad_test[2] * x56 + x52) +
                               grad_trial[1] * (grad_test[0] * x74 - grad_test[2] * (x66 + x75) + x71) +
                               grad_trial[2] * (grad_test[0] * x64 + grad_test[1] * x70 - x62));
                bf[4] += dx * (grad_trial[0] * (grad_test[1] * x43 + grad_test[2] * x49 - x25) -
                               grad_trial[1] * (grad_test[0] * x30 - grad_test[2] * x47 + x37) -
                               grad_trial[2] * (-grad_test[0] * x34 + grad_test[1] * x41 + x45));
                bf[5] += -x21 * (grad_trial[0] * (-grad_test[0] * (C1 + C2 * (x81 + x83)) + x20 * x77 + x6 * x79) +
                                 grad_trial[1] * (-grad_test[1] * (C1 + C2 * (x81 + x87)) + x19 * x77 + x6 * x85) +
                                 grad_trial[2] * (-grad_test[2] * (C1 + C2 * (x83 + x87)) + x19 * x79 + x20 * x85));
                bf[6] += dx * (-grad_trial[0] * (grad_test[1] * x94 - grad_test[2] * x98 + x90) -
                               grad_trial[1] * (-grad_test[0] * x106 + grad_test[2] * x104 + x100) +
                               grad_trial[2] * (grad_test[0] * x111 + grad_test[1] * x109 - x107));
                bf[8] += dx * (grad_trial[0] * (grad_test[1] * x74 + grad_test[2] * x64 - x52) +
                               grad_trial[1] * (-grad_test[0] * (x60 + x72) + grad_test[2] * x70 + x71) -
                               grad_trial[2] * (grad_test[0] * x56 - grad_test[1] * (x65 - x75) + x62));
                bf[9] += dx * (grad_trial[0] * (grad_test[1] * x106 + grad_test[2] * x111 - x90) -
                               grad_trial[1] * (grad_test[0] * x94 - grad_test[2] * x109 + x100) -
                               grad_trial[2] * (-grad_test[0] * x98 + grad_test[1] * x104 + x107));
                bf[10] +=
                    -x21 * (grad_trial[0] * (-grad_test[0] * (C1 + C2 * (x114 + x115)) + x112 * x20 + x113 * x6) +
                            grad_trial[1] * (-grad_test[1] * (C1 + C2 * (x114 + x117)) + x112 * x19 + x116 * x6) +
                            grad_trial[2] * (-grad_test[2] * (C1 + C2 * (x115 + x117)) + x113 * x19 + x116 * x20));
                bf[15] += 0;
                bf[12] += x121 * (grad_trial[0] * x118 - grad_trial[1] * x120 + grad_trial[2] * x119);
                bf[13] += -x121 * (grad_trial[0] * x122 - grad_trial[1] * x124 + grad_trial[2] * x123);
                bf[14] += x121 * (grad_trial[0] * x125 - grad_trial[1] * x127 + grad_trial[2] * x126);
                bf[3] += x128 * (grad_test[0] * x118 - grad_test[1] * x120 + grad_test[2] * x119);
                bf[7] += -x128 * (grad_test[0] * x122 - grad_test[1] * x124 + grad_test[2] * x123);
                bf[11] += x128 * (grad_test[0] * x125 - grad_test[1] * x127 + grad_test[2] * x126);
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T p,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T fun_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = 2 * C1;
                T x1 = f[4] * f[8];
                T x2 = f[5] * f[7];
                T x3 = pow(f[0], 2) + pow(f[3], 2) + pow(f[6], 2);
                T x4 = f[0] * f[1] + f[3] * f[4] + f[6] * f[7];
                T x5 = f[0] * f[2] + f[3] * f[5] + f[6] * f[8];
                T x6 = pow(f[2], 2) + pow(f[5], 2) + pow(f[8], 2);
                T x7 = pow(f[1], 2) + pow(f[4], 2) + pow(f[7], 2);
                T x8 = x3 + x6 + x7;
                T x9 = 2 * C2;
                T x10 = f[3] * f[7];
                T x11 = f[4] * f[6];
                T x12 = f[1] * f[2] + f[4] * f[5] + f[7] * f[8];
                T x13 = f[3] * f[8];
                T x14 = f[5] * f[6];
                lf[0] +=
                    dx *
                    (grad_test[0] * (f[0] * x0 + p * (x1 - x2) - x9 * (f[0] * x3 - f[0] * x8 + f[1] * x4 + f[2] * x5)) -
                     grad_test[1] *
                         (-f[1] * x0 + p * (x13 - x14) + x9 * (f[0] * x4 + f[1] * x7 - f[1] * x8 + f[2] * x12)) +
                     grad_test[2] *
                         (f[2] * x0 + p * (x10 - x11) - x9 * (f[0] * x5 + f[1] * x12 + f[2] * x6 - f[2] * x8)));
                lf[1] += dx * (-grad_test[0] * (-f[3] * x0 + p * (f[1] * f[8] - f[2] * f[7]) +
                                                x9 * (f[3] * x3 - f[3] * x8 + f[4] * x4 + f[5] * x5)) +
                               grad_test[1] * (f[4] * x0 + p * (f[0] * f[8] - f[2] * f[6]) -
                                               x9 * (f[3] * x4 + f[4] * x7 - f[4] * x8 + f[5] * x12)) -
                               grad_test[2] * (-f[5] * x0 + p * (f[0] * f[7] - f[1] * f[6]) +
                                               x9 * (f[3] * x5 + f[4] * x12 + f[5] * x6 - f[5] * x8)));
                lf[2] += dx * (grad_test[0] * (f[6] * x0 + p * (f[1] * f[5] - f[2] * f[4]) -
                                               x9 * (f[6] * x3 - f[6] * x8 + f[7] * x4 + f[8] * x5)) -
                               grad_test[1] * (-f[7] * x0 + p * (f[0] * f[5] - f[2] * f[3]) +
                                               x9 * (f[6] * x4 + f[7] * x7 - f[7] * x8 + f[8] * x12)) +
                               grad_test[2] * (f[8] * x0 + p * (f[0] * f[4] - f[1] * f[3]) -
                                               x9 * (f[6] * x5 + f[7] * x12 + f[8] * x6 - f[8] * x8)));
                lf[3] +=
                    dx * fun_test * (f[0] * x1 - f[0] * x2 - f[1] * x13 + f[1] * x14 + f[2] * x10 - f[2] * x11 - 1);
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
                T x3 = x1 + x2;
                T x4 = x0 + x3;
                T x5 = pow(f[1], 2);
                T x6 = pow(f[4], 2);
                T x7 = pow(f[7], 2);
                T x8 = x6 + x7;
                T x9 = x5 + x8;
                T x10 = pow(f[2], 2);
                T x11 = pow(f[5], 2);
                T x12 = pow(f[8], 2);
                T x13 = x11 + x12;
                T x14 = x10 + x13;
                T x15 = x14 + x4 + x9;
                T x16 = f[4] * f[8];
                T x17 = f[5] * f[6];
                T x18 = f[3] * f[7];
                T x19 = f[5] * f[7];
                T x20 = f[3] * f[8];
                T x21 = f[4] * f[6];
                T x22 = f[0] * x16 - f[0] * x19 + f[1] * x17 - f[1] * x20 + f[2] * x18 - f[2] * x21 - 1;
                T x23 = f[0] * f[1];
                T x24 = f[3] * f[4];
                T x25 = f[6] * f[7];
                T x26 = x24 + x25;
                T x27 = x23 + x26;
                T x28 = f[0] * f[2];
                T x29 = f[3] * f[5];
                T x30 = f[6] * f[8];
                T x31 = x29 + x30;
                T x32 = x28 + x31;
                T x33 = f[1] * f[2];
                T x34 = f[4] * f[5];
                T x35 = f[7] * f[8];
                T x36 = x34 + x35;
                T x37 = x33 + x36;
                T x38 = 2 * C1;
                T x39 = -x19;
                T x40 = x16 + x39;
                T x41 = 2 * C2;
                T x42 = -x21;
                T x43 = x18 + x42;
                T x44 = -x17;
                T x45 = x20 + x44;
                T x46 = f[0] * f[8];
                T x47 = f[2] * f[6];
                T x48 = -x47;
                T x49 = x46 + x48;
                T x50 = f[1] * f[8];
                T x51 = f[2] * f[7];
                T x52 = -x51;
                T x53 = x50 + x52;
                T x54 = f[0] * f[7];
                T x55 = f[1] * f[6];
                T x56 = -x55;
                T x57 = x54 + x56;
                T x58 = f[1] * f[5];
                T x59 = f[2] * f[4];
                T x60 = -x59;
                T x61 = x58 + x60;
                T x62 = f[0] * f[4];
                T x63 = f[1] * f[3];
                T x64 = -x63;
                T x65 = x62 + x64;
                T x66 = f[0] * f[5];
                T x67 = f[2] * f[3];
                T x68 = -x67;
                T x69 = x66 + x68;
                T x70 = dx * fun_test;
                T x71 = C2 * x26;
                T x72 = C2 * grad_test[2];
                T x73 = C2 * grad_test[0];
                T x74 = C2 * grad_test[1];
                T x75 = 2 * dx;
                T x76 = f[1] * f[4];
                T x77 = f[2] * f[5];
                T x78 = grad_test[0] * x41;
                T x79 = x78 * (x76 + x77);
                T x80 = f[8] * p;
                T x81 = x41 * (x62 - 2 * x63) + x80;
                T x82 = f[7] * p;
                T x83 = -x41 * (x66 - 2 * x67) + x82;
                T x84 = f[0] * f[3];
                T x85 = grad_test[1] * x41;
                T x86 = x85 * (x77 + x84);
                T x87 = f[6] * p;
                T x88 = x41 * (x58 - 2 * x59) + x87;
                T x89 = x41 * (2 * x62 + x64) + x80;
                T x90 = grad_test[2] * x41;
                T x91 = x90 * (x76 + x84);
                T x92 = x41 * (2 * x58 + x60) + x87;
                T x93 = x41 * (2 * x66 + x68) - x82;
                T x94 = f[1] * f[7];
                T x95 = f[2] * f[8];
                T x96 = x78 * (x94 + x95);
                T x97 = f[4] * p;
                T x98 = x41 * (x46 - 2 * x47) + x97;
                T x99 = f[5] * p;
                T x100 = x41 * (x54 - 2 * x55);
                T x101 = f[0] * f[6];
                T x102 = x90 * (x101 + x94);
                T x103 = x41 * (2 * x46 + x48) + x97;
                T x104 = f[3] * p;
                T x105 = -x104;
                T x106 = x105 + x41 * (2 * x50 + x52);
                T x107 = -x85 * (x101 + x95);
                T x108 = -x99;
                T x109 = x108 + x41 * (2 * x54 + x56);
                T x110 = x41 * (x50 - 2 * x51);
                T x111 = x23 + x25;
                T x112 = x28 + x30;
                T x113 = x10 + x12;
                T x114 = x5 + x7;
                T x115 = x33 + x35;
                T x116 = x0 + x2;
                T x117 = f[4] * f[7];
                T x118 = f[5] * f[8];
                T x119 = x78 * (x117 + x118);
                T x120 = f[2] * p;
                T x121 = x120 + x41 * (x18 - 2 * x21);
                T x122 = f[1] * p;
                T x123 = x122 - x41 * (-2 * x17 + x20);
                T x124 = f[3] * f[6];
                T x125 = x85 * (x118 + x124);
                T x126 = f[0] * p;
                T x127 = x126 + x41 * (x16 - 2 * x19);
                T x128 = x120 + x41 * (2 * x18 + x42);
                T x129 = x90 * (x117 + x124);
                T x130 = x126 + x41 * (2 * x16 + x39);
                T x131 = -x122 + x41 * (2 * x20 + x44);
                T x132 = x23 + x24;
                T x133 = x28 + x29;
                T x134 = x10 + x11;
                T x135 = x5 + x6;
                T x136 = x33 + x34;
                T x137 = x0 + x1;
                T x138 = dx * fun_trial;
                e += dx * (C1 * (x15 - 3) +
                           C2 * (-1.0 / 2.0 * pow(x14, 2) + (1.0 / 2.0) * pow(x15, 2) - pow(x27, 2) - pow(x32, 2) -
                                 pow(x37, 2) - 1.0 / 2.0 * pow(x4, 2) - 1.0 / 2.0 * pow(x9, 2) - 3) +
                           p * x22);
                lf[0] +=
                    dx *
                    (grad_test[0] * (f[0] * x38 + p * x40 - x41 * (-f[0] * x15 + f[0] * x4 + f[1] * x27 + f[2] * x32)) -
                     grad_test[1] * (-f[1] * x38 + p * x45 + x41 * (f[0] * x27 - f[1] * x15 + f[1] * x9 + f[2] * x37)) +
                     grad_test[2] * (f[2] * x38 + p * x43 - x41 * (f[0] * x32 + f[1] * x37 + f[2] * x14 - f[2] * x15)));
                lf[1] += dx * (-grad_test[0] *
                                   (-f[3] * x38 + p * x53 + x41 * (-f[3] * x15 + f[3] * x4 + f[4] * x27 + f[5] * x32)) +
                               grad_test[1] *
                                   (f[4] * x38 + p * x49 - x41 * (f[3] * x27 - f[4] * x15 + f[4] * x9 + f[5] * x37)) -
                               grad_test[2] *
                                   (-f[5] * x38 + p * x57 + x41 * (f[3] * x32 + f[4] * x37 + f[5] * x14 - f[5] * x15)));
                lf[2] +=
                    dx *
                    (grad_test[0] * (f[6] * x38 + p * x61 - x41 * (-f[6] * x15 + f[6] * x4 + f[7] * x27 + f[8] * x32)) -
                     grad_test[1] * (-f[7] * x38 + p * x69 + x41 * (f[6] * x27 - f[7] * x15 + f[7] * x9 + f[8] * x37)) +
                     grad_test[2] * (f[8] * x38 + p * x65 - x41 * (f[6] * x32 + f[7] * x37 + f[8] * x14 - f[8] * x15)));
                lf[3] += x22 * x70;
                bf[0] +=
                    -x75 * (grad_trial[0] * (-grad_test[0] * (C1 + C2 * (x13 + x8)) + grad_test[1] * x71 + x31 * x72) +
                            grad_trial[1] * (grad_test[0] * x71 - grad_test[1] * (C1 + C2 * (x13 + x3)) + x36 * x72) +
                            grad_trial[2] * (-grad_test[2] * (C1 + C2 * (x3 + x8)) + x31 * x73 + x36 * x74));
                bf[1] += dx * (-grad_trial[0] * (grad_test[1] * x81 - grad_test[2] * x83 + x79) -
                               grad_trial[1] * (-grad_test[0] * x89 + grad_test[2] * x88 + x86) +
                               grad_trial[2] * (grad_test[0] * x93 + grad_test[1] * x92 - x91));
                bf[2] += dx * (-grad_trial[0] * (-grad_test[1] * (-x100 + x99) + grad_test[2] * x98 + x96) +
                               grad_trial[1] * (grad_test[0] * x109 - grad_test[2] * (x105 + x110) + x107) +
                               grad_trial[2] * (grad_test[0] * x103 + grad_test[1] * x106 - x102));
                bf[4] += dx * (grad_trial[0] * (grad_test[1] * x89 + grad_test[2] * x93 - x79) -
                               grad_trial[1] * (grad_test[0] * x81 - grad_test[2] * x92 + x86) -
                               grad_trial[2] * (-grad_test[0] * x83 + grad_test[1] * x88 + x91));
                bf[5] += -x75 * (grad_trial[0] * (-grad_test[0] * (C1 + C2 * (x113 + x114)) + x111 * x74 + x112 * x72) +
                                 grad_trial[1] * (-grad_test[1] * (C1 + C2 * (x113 + x116)) + x111 * x73 + x115 * x72) +
                                 grad_trial[2] * (-grad_test[2] * (C1 + C2 * (x114 + x116)) + x112 * x73 + x115 * x74));
                bf[6] += dx * (-grad_trial[0] * (grad_test[1] * x121 - grad_test[2] * x123 + x119) -
                               grad_trial[1] * (-grad_test[0] * x128 + grad_test[2] * x127 + x125) +
                               grad_trial[2] * (grad_test[0] * x131 + grad_test[1] * x130 - x129));
                bf[8] += dx * (grad_trial[0] * (grad_test[1] * x109 + grad_test[2] * x103 - x96) +
                               grad_trial[1] * (-grad_test[0] * (x100 + x108) + grad_test[2] * x106 + x107) -
                               grad_trial[2] * (grad_test[0] * x98 - grad_test[1] * (x104 - x110) + x102));
                bf[9] += dx * (grad_trial[0] * (grad_test[1] * x128 + grad_test[2] * x131 - x119) -
                               grad_trial[1] * (grad_test[0] * x121 - grad_test[2] * x130 + x125) -
                               grad_trial[2] * (-grad_test[0] * x123 + grad_test[1] * x127 + x129));
                bf[10] +=
                    -x75 * (grad_trial[0] * (-grad_test[0] * (C1 + C2 * (x134 + x135)) + x132 * x74 + x133 * x72) +
                            grad_trial[1] * (-grad_test[1] * (C1 + C2 * (x134 + x137)) + x132 * x73 + x136 * x72) +
                            grad_trial[2] * (-grad_test[2] * (C1 + C2 * (x135 + x137)) + x133 * x73 + x136 * x74));
                bf[15] += 0;
                bf[12] += x70 * (grad_trial[0] * x40 - grad_trial[1] * x45 + grad_trial[2] * x43);
                bf[13] += -x70 * (grad_trial[0] * x53 - grad_trial[1] * x49 + grad_trial[2] * x57);
                bf[14] += x70 * (grad_trial[0] * x61 - grad_trial[1] * x69 + grad_trial[2] * x65);
                bf[3] += x138 * (grad_test[0] * x40 - grad_test[1] * x45 + grad_test[2] * x43);
                bf[7] += -x138 * (grad_test[0] * x53 - grad_test[1] * x49 + grad_test[2] * x57);
                bf[11] += x138 * (grad_test[0] * x61 - grad_test[1] * x69 + grad_test[2] * x65);
            }

            T C1{1.0};
            T C2{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_3_IMPL_hpp
