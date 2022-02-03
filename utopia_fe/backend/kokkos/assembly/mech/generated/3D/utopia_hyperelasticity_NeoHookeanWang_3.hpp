#ifndef UTOPIA_TPL_HYPERELASTICITY_NeoHookeanWang_3_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeoHookeanWang_3_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_NeoHookeanWang.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of NeoHookeanWang for dimension 3
         */
        template <typename T>
        class NeoHookeanWang<T, 3> {
        public:
            static constexpr int Dim = 3;

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    StressStrainParameters<T, T> ssp;
                    ssp.read(in);

                    lambda = ssp.first_lame_parameter.get();
                    mu = ssp.shear_modulus.get();
                }

                T mu{1.0};
                T lambda{1.0};
            };

            NeoHookeanWang(const Params &params) {
                mu = params.mu;
                lambda = params.lambda;
            }

            UTOPIA_FUNCTION void hessian(const T *UTOPIA_RESTRICT f,
                                         const T *grad_test,
                                         const T *grad_trial,
                                         const T dx,
                                         T *UTOPIA_RESTRICT bf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[4] * f[8];
                T x1 = f[5] * f[6];
                T x2 = f[3] * f[7];
                T x3 = f[5] * f[7];
                T x4 = f[3] * f[8];
                T x5 = f[4] * f[6];
                T x6 = f[0] * x0 - f[0] * x3 + f[1] * x1 - f[1] * x4 + f[2] * x2 - f[2] * x5;
                T x7 = 2 / pow(x6, 2.0 / 3.0);
                T x8 = -2.0 / 3.0 * x0 + (2.0 / 3.0) * x3;
                T x9 = pow(x6, -5.0 / 3.0);
                T x10 = f[0] * x9;
                T x11 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) +
                        pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2);
                T x12 = x11 / pow(x6, 8.0 / 3.0);
                T x13 = x12 * (-5.0 / 3.0 * x0 + (5.0 / 3.0) * x3);
                T x14 = (1.0 / 2.0) * mu;
                T x15 = grad_test[0] * x14;
                T x16 = -2.0 / 3.0 * x1 + (2.0 / 3.0) * x4;
                T x17 = 2 * x10;
                T x18 = f[1] * x9;
                T x19 = 2 * x8;
                T x20 = x16 * x17 + x18 * x19;
                T x21 = grad_test[1] * x14;
                T x22 = -2.0 / 3.0 * x2 + (2.0 / 3.0) * x5;
                T x23 = f[2] * x9;
                T x24 = x17 * x22 + x19 * x23;
                T x25 = grad_test[2] * x14;
                T x26 = x12 * (-5.0 / 3.0 * x1 + (5.0 / 3.0) * x4);
                T x27 = 2 * x18;
                T x28 = 2 * x16;
                T x29 = x22 * x27 + x23 * x28;
                T x30 = x12 * (-5.0 / 3.0 * x2 + (5.0 / 3.0) * x5);
                T x31 = f[1] * f[8];
                T x32 = f[2] * f[7];
                T x33 = x12 * ((5.0 / 3.0) * x31 - 5.0 / 3.0 * x32);
                T x34 = (2.0 / 3.0) * x31 - 2.0 / 3.0 * x32;
                T x35 = f[3] * x9;
                T x36 = x17 * x34 + x19 * x35;
                T x37 = (1.0 / 2.0) * lambda;
                T x38 = f[8] * x37;
                T x39 = -x38;
                T x40 = (2.0 / 3.0) * x11;
                T x41 = f[8] * x9;
                T x42 = x40 * x41;
                T x43 = x27 * x34 + x28 * x35 + x42;
                T x44 = f[7] * x37;
                T x45 = 2 * x23;
                T x46 = 2 * x22;
                T x47 = x40 * x9;
                T x48 = f[7] * x47;
                T x49 = x34 * x45 + x35 * x46 - x48;
                T x50 = f[0] * f[8];
                T x51 = f[2] * f[6];
                T x52 = x12 * (-5.0 / 3.0 * x50 + (5.0 / 3.0) * x51);
                T x53 = -2.0 / 3.0 * x50 + (2.0 / 3.0) * x51;
                T x54 = f[4] * x9;
                T x55 = x27 * x53 + x28 * x54;
                T x56 = x17 * x53 + x19 * x54 - x42;
                T x57 = f[6] * x37;
                T x58 = -x57;
                T x59 = f[6] * x47;
                T x60 = x45 * x53 + x46 * x54 + x59;
                T x61 = f[0] * f[7];
                T x62 = f[1] * f[6];
                T x63 = x12 * ((5.0 / 3.0) * x61 - 5.0 / 3.0 * x62);
                T x64 = (2.0 / 3.0) * x61 - 2.0 / 3.0 * x62;
                T x65 = f[5] * x9;
                T x66 = x45 * x64 + x46 * x65;
                T x67 = -x44;
                T x68 = x17 * x64 + x19 * x65 + x48;
                T x69 = x27 * x64 + x28 * x65 - x59;
                T x70 = f[1] * f[5];
                T x71 = f[2] * f[4];
                T x72 = x12 * (-5.0 / 3.0 * x70 + (5.0 / 3.0) * x71);
                T x73 = -2.0 / 3.0 * x70 + (2.0 / 3.0) * x71;
                T x74 = f[6] * x9;
                T x75 = x17 * x73 + x19 * x74;
                T x76 = f[5] * x37;
                T x77 = f[5] * x47;
                T x78 = x27 * x73 + x28 * x74 - x77;
                T x79 = f[4] * x37;
                T x80 = -x79;
                T x81 = f[4] * x47;
                T x82 = x45 * x73 + x46 * x74 + x81;
                T x83 = f[0] * f[5];
                T x84 = f[2] * f[3];
                T x85 = x12 * ((5.0 / 3.0) * x83 - 5.0 / 3.0 * x84);
                T x86 = (2.0 / 3.0) * x83 - 2.0 / 3.0 * x84;
                T x87 = f[7] * x9;
                T x88 = x27 * x86 + x28 * x87;
                T x89 = -x76;
                T x90 = x17 * x86 + x19 * x87 + x77;
                T x91 = f[3] * x37;
                T x92 = x35 * x40;
                T x93 = x45 * x86 + x46 * x87 - x92;
                T x94 = f[0] * f[4];
                T x95 = f[1] * f[3];
                T x96 = x12 * (-5.0 / 3.0 * x94 + (5.0 / 3.0) * x95);
                T x97 = -2.0 / 3.0 * x94 + (2.0 / 3.0) * x95;
                T x98 = x41 * x46 + x45 * x97;
                T x99 = x17 * x97 + x19 * x41 - x81;
                T x100 = -x91;
                T x101 = x27 * x97 + x28 * x41 + x92;
                T x102 = 2 * x35;
                T x103 = 2 * x34;
                T x104 = x102 * x53 + x103 * x54;
                T x105 = x102 * x64 + x103 * x65;
                T x106 = 2 * x54;
                T x107 = 2 * x53;
                T x108 = x106 * x64 + x107 * x65;
                T x109 = x102 * x73 + x103 * x74;
                T x110 = f[2] * x37;
                T x111 = -x110;
                T x112 = x23 * x40;
                T x113 = x106 * x73 + x107 * x74 + x112;
                T x114 = f[1] * x37;
                T x115 = 2 * x65;
                T x116 = 2 * x64;
                T x117 = x18 * x40;
                T x118 = x115 * x73 + x116 * x74 - x117;
                T x119 = x106 * x86 + x107 * x87;
                T x120 = x102 * x86 + x103 * x87 - x112;
                T x121 = f[0] * x37;
                T x122 = -x121;
                T x123 = x10 * x40;
                T x124 = x115 * x86 + x116 * x87 + x123;
                T x125 = x115 * x97 + x116 * x41;
                T x126 = -x114;
                T x127 = x102 * x97 + x103 * x41 + x117;
                T x128 = x106 * x97 + x107 * x41 - x123;
                T x129 = 2 * x74;
                T x130 = 2 * x73;
                T x131 = x129 * x86 + x130 * x87;
                T x132 = x129 * x97 + x130 * x41;
                T x133 = 2 * x41 * x86 + 2 * x87 * x97;
                bf[0] += dx * (grad_trial[0] * (x15 * (4 * x10 * x8 + x13 * x8 + x7) + x21 * (x13 * x16 + x20) +
                                                x25 * (x13 * x22 + x24)) +
                               grad_trial[1] * (x15 * (x20 + x26 * x8) + x21 * (4 * x16 * x18 + x16 * x26 + x7) +
                                                x25 * (x22 * x26 + x29)) +
                               grad_trial[2] * (x15 * (x24 + x30 * x8) + x21 * (x16 * x30 + x29) +
                                                x25 * (4 * x22 * x23 + x22 * x30 + x7)));
                bf[1] +=
                    dx * (grad_trial[0] * (grad_test[1] * (x14 * (x16 * x33 + x43) + x39) +
                                           grad_test[2] * (x14 * (x22 * x33 + x49) + x44) + x15 * (x33 * x8 + x36)) +
                          grad_trial[1] * (grad_test[0] * (x14 * (x52 * x8 + x56) + x38) +
                                           grad_test[2] * (x14 * (x22 * x52 + x60) + x58) + x21 * (x16 * x52 + x55)) +
                          grad_trial[2] * (grad_test[0] * (x14 * (x63 * x8 + x68) + x67) +
                                           grad_test[1] * (x14 * (x16 * x63 + x69) + x57) + x25 * (x22 * x63 + x66)));
                bf[2] +=
                    dx * (grad_trial[0] * (grad_test[1] * (x14 * (x16 * x72 + x78) + x76) +
                                           grad_test[2] * (x14 * (x22 * x72 + x82) + x80) + x15 * (x72 * x8 + x75)) +
                          grad_trial[1] * (grad_test[0] * (x14 * (x8 * x85 + x90) + x89) +
                                           grad_test[2] * (x14 * (x22 * x85 + x93) + x91) + x21 * (x16 * x85 + x88)) +
                          grad_trial[2] * (grad_test[0] * (x14 * (x8 * x96 + x99) + x79) +
                                           grad_test[1] * (x100 + x14 * (x101 + x16 * x96)) + x25 * (x22 * x96 + x98)));
                bf[3] +=
                    dx * (grad_trial[0] * (grad_test[1] * (x14 * (x13 * x53 + x56) + x38) +
                                           grad_test[2] * (x14 * (x13 * x64 + x68) + x67) + x15 * (x13 * x34 + x36)) +
                          grad_trial[1] * (grad_test[0] * (x14 * (x26 * x34 + x43) + x39) +
                                           grad_test[2] * (x14 * (x26 * x64 + x69) + x57) + x21 * (x26 * x53 + x55)) +
                          grad_trial[2] * (grad_test[0] * (x14 * (x30 * x34 + x49) + x44) +
                                           grad_test[1] * (x14 * (x30 * x53 + x60) + x58) + x25 * (x30 * x64 + x66)));
                bf[4] += dx * (grad_trial[0] * (x15 * (x33 * x34 + 4 * x34 * x35 + x7) + x21 * (x104 + x33 * x53) +
                                                x25 * (x105 + x33 * x64)) +
                               grad_trial[1] * (x15 * (x104 + x34 * x52) + x21 * (x52 * x53 + 4 * x53 * x54 + x7) +
                                                x25 * (x108 + x52 * x64)) +
                               grad_trial[2] * (x15 * (x105 + x34 * x63) + x21 * (x108 + x53 * x63) +
                                                x25 * (x63 * x64 + 4 * x64 * x65 + x7)));
                bf[5] +=
                    dx *
                    (grad_trial[0] * (grad_test[1] * (x111 + x14 * (x113 + x53 * x72)) +
                                      grad_test[2] * (x114 + x14 * (x118 + x64 * x72)) + x15 * (x109 + x34 * x72)) +
                     grad_trial[1] * (grad_test[0] * (x110 + x14 * (x120 + x34 * x85)) +
                                      grad_test[2] * (x122 + x14 * (x124 + x64 * x85)) + x21 * (x119 + x53 * x85)) +
                     grad_trial[2] * (grad_test[0] * (x126 + x14 * (x127 + x34 * x96)) +
                                      grad_test[1] * (x121 + x14 * (x128 + x53 * x96)) + x25 * (x125 + x64 * x96)));
                bf[6] +=
                    dx * (grad_trial[0] * (grad_test[1] * (x14 * (x13 * x86 + x90) + x89) +
                                           grad_test[2] * (x14 * (x13 * x97 + x99) + x79) + x15 * (x13 * x73 + x75)) +
                          grad_trial[1] * (grad_test[0] * (x14 * (x26 * x73 + x78) + x76) +
                                           grad_test[2] * (x100 + x14 * (x101 + x26 * x97)) + x21 * (x26 * x86 + x88)) +
                          grad_trial[2] * (grad_test[0] * (x14 * (x30 * x73 + x82) + x80) +
                                           grad_test[1] * (x14 * (x30 * x86 + x93) + x91) + x25 * (x30 * x97 + x98)));
                bf[7] +=
                    dx *
                    (grad_trial[0] * (grad_test[1] * (x110 + x14 * (x120 + x33 * x86)) +
                                      grad_test[2] * (x126 + x14 * (x127 + x33 * x97)) + x15 * (x109 + x33 * x73)) +
                     grad_trial[1] * (grad_test[0] * (x111 + x14 * (x113 + x52 * x73)) +
                                      grad_test[2] * (x121 + x14 * (x128 + x52 * x97)) + x21 * (x119 + x52 * x86)) +
                     grad_trial[2] * (grad_test[0] * (x114 + x14 * (x118 + x63 * x73)) +
                                      grad_test[1] * (x122 + x14 * (x124 + x63 * x86)) + x25 * (x125 + x63 * x97)));
                bf[8] += dx * (grad_trial[0] * (x15 * (x7 + x72 * x73 + 4 * x73 * x74) + x21 * (x131 + x72 * x86) +
                                                x25 * (x132 + x72 * x97)) +
                               grad_trial[1] * (x15 * (x131 + x73 * x85) + x21 * (x7 + x85 * x86 + 4 * x86 * x87) +
                                                x25 * (x133 + x85 * x97)) +
                               grad_trial[2] * (x15 * (x132 + x73 * x96) + x21 * (x133 + x86 * x96) +
                                                x25 * (4 * x41 * x97 + x7 + x96 * x97)));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[4] * f[8];
                T x1 = f[5] * f[7];
                T x2 = (1.0 / 2.0) * lambda;
                T x3 = f[5] * f[6];
                T x4 = f[3] * f[7];
                T x5 = f[3] * f[8];
                T x6 = f[4] * f[6];
                T x7 = f[0] * x0 - f[0] * x1 + f[1] * x3 - f[1] * x5 + f[2] * x4 - f[2] * x6;
                T x8 = 2 / pow(x7, 2.0 / 3.0);
                T x9 = (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) +
                        pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2)) /
                       pow(x7, 5.0 / 3.0);
                T x10 = (1.0 / 2.0) * mu;
                T x11 = f[2] * f[7];
                T x12 = f[1] * f[8];
                T x13 = f[0] * f[8];
                T x14 = f[2] * f[6];
                T x15 = f[1] * f[6];
                T x16 = f[0] * f[7];
                T x17 = f[1] * f[5];
                T x18 = f[2] * f[4];
                T x19 = f[2] * f[3];
                T x20 = f[0] * f[5];
                T x21 = f[0] * f[4];
                T x22 = f[1] * f[3];
                lf[0] +=
                    dx *
                    (grad_test[0] * (x10 * (f[0] * x8 + x9 * (-2.0 / 3.0 * x0 + (2.0 / 3.0) * x1)) + x2 * (x0 - x1)) +
                     grad_test[1] * (x10 * (f[1] * x8 + x9 * (-2.0 / 3.0 * x3 + (2.0 / 3.0) * x5)) + x2 * (x3 - x5)) +
                     grad_test[2] * (x10 * (f[2] * x8 + x9 * (-2.0 / 3.0 * x4 + (2.0 / 3.0) * x6)) + x2 * (x4 - x6)));
                lf[1] +=
                    dx * (grad_test[0] *
                              (x10 * (f[3] * x8 + x9 * (-2.0 / 3.0 * x11 + (2.0 / 3.0) * x12)) + x2 * (x11 - x12)) +
                          grad_test[1] *
                              (x10 * (f[4] * x8 + x9 * (-2.0 / 3.0 * x13 + (2.0 / 3.0) * x14)) + x2 * (x13 - x14)) +
                          grad_test[2] *
                              (x10 * (f[5] * x8 + x9 * (-2.0 / 3.0 * x15 + (2.0 / 3.0) * x16)) + x2 * (x15 - x16)));
                lf[2] +=
                    dx * (grad_test[0] *
                              (x10 * (f[6] * x8 + x9 * (-2.0 / 3.0 * x17 + (2.0 / 3.0) * x18)) + x2 * (x17 - x18)) +
                          grad_test[1] *
                              (x10 * (f[7] * x8 + x9 * (-2.0 / 3.0 * x19 + (2.0 / 3.0) * x20)) + x2 * (x19 - x20)) +
                          grad_test[2] *
                              (x10 * (f[8] * x8 + x9 * (-2.0 / 3.0 * x21 + (2.0 / 3.0) * x22)) + x2 * (x21 - x22)));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[0] * f[4] * f[8] - f[0] * f[5] * f[7] - f[1] * f[3] * f[8] + f[1] * f[5] * f[6] +
                       f[2] * f[3] * f[7] - f[2] * f[4] * f[6];
                e += dx * ((1.0 / 2.0) * lambda * (x0 - 1) +
                           (1.0 / 2.0) * mu *
                               (-3 + (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) +
                                      pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2)) /
                                         pow(x0, 2.0 / 3.0)));
            }

            UTOPIA_FUNCTION void eval(const T *UTOPIA_RESTRICT f,
                                      const T *grad_test,
                                      const T *grad_trial,
                                      const T dx,
                                      T &e,
                                      T *UTOPIA_RESTRICT lf,
                                      T *UTOPIA_RESTRICT bf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[4] * f[8];
                T x1 = f[5] * f[6];
                T x2 = f[3] * f[7];
                T x3 = f[5] * f[7];
                T x4 = f[3] * f[8];
                T x5 = f[4] * f[6];
                T x6 = f[0] * x0 - f[0] * x3 + f[1] * x1 - f[1] * x4 + f[2] * x2 - f[2] * x5;
                T x7 = (1.0 / 2.0) * lambda;
                T x8 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) +
                       pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2);
                T x9 = pow(x6, -2.0 / 3.0);
                T x10 = (1.0 / 2.0) * mu;
                T x11 = 2 * x9;
                T x12 = -2.0 / 3.0 * x0 + (2.0 / 3.0) * x3;
                T x13 = pow(x6, -5.0 / 3.0);
                T x14 = x13 * x8;
                T x15 = -2.0 / 3.0 * x1 + (2.0 / 3.0) * x4;
                T x16 = -2.0 / 3.0 * x2 + (2.0 / 3.0) * x5;
                T x17 = f[2] * f[7];
                T x18 = f[1] * f[8];
                T x19 = -2.0 / 3.0 * x17 + (2.0 / 3.0) * x18;
                T x20 = f[0] * f[8];
                T x21 = f[2] * f[6];
                T x22 = -2.0 / 3.0 * x20 + (2.0 / 3.0) * x21;
                T x23 = f[1] * f[6];
                T x24 = f[0] * f[7];
                T x25 = -2.0 / 3.0 * x23 + (2.0 / 3.0) * x24;
                T x26 = f[1] * f[5];
                T x27 = f[2] * f[4];
                T x28 = -2.0 / 3.0 * x26 + (2.0 / 3.0) * x27;
                T x29 = f[2] * f[3];
                T x30 = f[0] * f[5];
                T x31 = -2.0 / 3.0 * x29 + (2.0 / 3.0) * x30;
                T x32 = f[0] * f[4];
                T x33 = f[1] * f[3];
                T x34 = -2.0 / 3.0 * x32 + (2.0 / 3.0) * x33;
                T x35 = f[0] * x13;
                T x36 = x8 / pow(x6, 8.0 / 3.0);
                T x37 = x36 * (-5.0 / 3.0 * x0 + (5.0 / 3.0) * x3);
                T x38 = grad_test[0] * x10;
                T x39 = 2 * x35;
                T x40 = f[1] * x13;
                T x41 = 2 * x12;
                T x42 = x15 * x39 + x40 * x41;
                T x43 = grad_test[1] * x10;
                T x44 = f[2] * x13;
                T x45 = x16 * x39 + x41 * x44;
                T x46 = grad_test[2] * x10;
                T x47 = x36 * (-5.0 / 3.0 * x1 + (5.0 / 3.0) * x4);
                T x48 = 2 * x40;
                T x49 = 2 * x15;
                T x50 = x16 * x48 + x44 * x49;
                T x51 = x36 * (-5.0 / 3.0 * x2 + (5.0 / 3.0) * x5);
                T x52 = x36 * (-5.0 / 3.0 * x17 + (5.0 / 3.0) * x18);
                T x53 = f[3] * x13;
                T x54 = x19 * x39 + x41 * x53;
                T x55 = f[8] * x7;
                T x56 = -x55;
                T x57 = (2.0 / 3.0) * x14;
                T x58 = f[8] * x57;
                T x59 = x19 * x48 + x49 * x53 + x58;
                T x60 = f[7] * x7;
                T x61 = 2 * x44;
                T x62 = 2 * x16;
                T x63 = f[7] * x57;
                T x64 = x19 * x61 + x53 * x62 - x63;
                T x65 = x36 * (-5.0 / 3.0 * x20 + (5.0 / 3.0) * x21);
                T x66 = f[4] * x13;
                T x67 = x22 * x48 + x49 * x66;
                T x68 = x22 * x39 + x41 * x66 - x58;
                T x69 = f[6] * x7;
                T x70 = -x69;
                T x71 = f[6] * x57;
                T x72 = x22 * x61 + x62 * x66 + x71;
                T x73 = x36 * (-5.0 / 3.0 * x23 + (5.0 / 3.0) * x24);
                T x74 = f[5] * x13;
                T x75 = x25 * x61 + x62 * x74;
                T x76 = -x60;
                T x77 = x25 * x39 + x41 * x74 + x63;
                T x78 = x25 * x48 + x49 * x74 - x71;
                T x79 = x36 * (-5.0 / 3.0 * x26 + (5.0 / 3.0) * x27);
                T x80 = f[6] * x13;
                T x81 = x28 * x39 + x41 * x80;
                T x82 = f[5] * x7;
                T x83 = f[5] * x57;
                T x84 = x28 * x48 + x49 * x80 - x83;
                T x85 = f[4] * x7;
                T x86 = -x85;
                T x87 = f[4] * x57;
                T x88 = x28 * x61 + x62 * x80 + x87;
                T x89 = x36 * (-5.0 / 3.0 * x29 + (5.0 / 3.0) * x30);
                T x90 = f[7] * x13;
                T x91 = x31 * x48 + x49 * x90;
                T x92 = -x82;
                T x93 = x31 * x39 + x41 * x90 + x83;
                T x94 = f[3] * x7;
                T x95 = f[3] * x57;
                T x96 = x31 * x61 + x62 * x90 - x95;
                T x97 = x36 * (-5.0 / 3.0 * x32 + (5.0 / 3.0) * x33);
                T x98 = f[8] * x13;
                T x99 = x34 * x61 + x62 * x98;
                T x100 = x34 * x39 + x41 * x98 - x87;
                T x101 = -x94;
                T x102 = x34 * x48 + x49 * x98 + x95;
                T x103 = 2 * x53;
                T x104 = 2 * x19;
                T x105 = x103 * x22 + x104 * x66;
                T x106 = x103 * x25 + x104 * x74;
                T x107 = 2 * x66;
                T x108 = 2 * x22;
                T x109 = x107 * x25 + x108 * x74;
                T x110 = x103 * x28 + x104 * x80;
                T x111 = f[2] * x7;
                T x112 = -x111;
                T x113 = f[2] * x57;
                T x114 = x107 * x28 + x108 * x80 + x113;
                T x115 = f[1] * x7;
                T x116 = 2 * x74;
                T x117 = 2 * x25;
                T x118 = f[1] * x57;
                T x119 = x116 * x28 + x117 * x80 - x118;
                T x120 = x107 * x31 + x108 * x90;
                T x121 = x103 * x31 + x104 * x90 - x113;
                T x122 = f[0] * x7;
                T x123 = -x122;
                T x124 = f[0] * x57;
                T x125 = x116 * x31 + x117 * x90 + x124;
                T x126 = x116 * x34 + x117 * x98;
                T x127 = -x115;
                T x128 = x103 * x34 + x104 * x98 + x118;
                T x129 = x107 * x34 + x108 * x98 - x124;
                T x130 = 2 * x80;
                T x131 = 2 * x28;
                T x132 = x130 * x31 + x131 * x90;
                T x133 = x130 * x34 + x131 * x98;
                T x134 = 2 * x31 * x98 + 2 * x34 * x90;
                e += dx * (x10 * (x8 * x9 - 3) + x7 * (x6 - 1));
                lf[0] += dx * (grad_test[0] * (x10 * (f[0] * x11 + x12 * x14) + x7 * (x0 - x3)) +
                               grad_test[1] * (x10 * (f[1] * x11 + x14 * x15) + x7 * (x1 - x4)) +
                               grad_test[2] * (x10 * (f[2] * x11 + x14 * x16) + x7 * (x2 - x5)));
                lf[1] += dx * (grad_test[0] * (x10 * (f[3] * x11 + x14 * x19) + x7 * (x17 - x18)) +
                               grad_test[1] * (x10 * (f[4] * x11 + x14 * x22) + x7 * (x20 - x21)) +
                               grad_test[2] * (x10 * (f[5] * x11 + x14 * x25) + x7 * (x23 - x24)));
                lf[2] += dx * (grad_test[0] * (x10 * (f[6] * x11 + x14 * x28) + x7 * (x26 - x27)) +
                               grad_test[1] * (x10 * (f[7] * x11 + x14 * x31) + x7 * (x29 - x30)) +
                               grad_test[2] * (x10 * (f[8] * x11 + x14 * x34) + x7 * (x32 - x33)));
                bf[0] += dx * (grad_trial[0] * (x38 * (x11 + 4 * x12 * x35 + x12 * x37) + x43 * (x15 * x37 + x42) +
                                                x46 * (x16 * x37 + x45)) +
                               grad_trial[1] * (x38 * (x12 * x47 + x42) + x43 * (x11 + 4 * x15 * x40 + x15 * x47) +
                                                x46 * (x16 * x47 + x50)) +
                               grad_trial[2] * (x38 * (x12 * x51 + x45) + x43 * (x15 * x51 + x50) +
                                                x46 * (x11 + 4 * x16 * x44 + x16 * x51)));
                bf[1] +=
                    dx * (grad_trial[0] * (grad_test[1] * (x10 * (x15 * x52 + x59) + x56) +
                                           grad_test[2] * (x10 * (x16 * x52 + x64) + x60) + x38 * (x12 * x52 + x54)) +
                          grad_trial[1] * (grad_test[0] * (x10 * (x12 * x65 + x68) + x55) +
                                           grad_test[2] * (x10 * (x16 * x65 + x72) + x70) + x43 * (x15 * x65 + x67)) +
                          grad_trial[2] * (grad_test[0] * (x10 * (x12 * x73 + x77) + x76) +
                                           grad_test[1] * (x10 * (x15 * x73 + x78) + x69) + x46 * (x16 * x73 + x75)));
                bf[2] +=
                    dx * (grad_trial[0] * (grad_test[1] * (x10 * (x15 * x79 + x84) + x82) +
                                           grad_test[2] * (x10 * (x16 * x79 + x88) + x86) + x38 * (x12 * x79 + x81)) +
                          grad_trial[1] * (grad_test[0] * (x10 * (x12 * x89 + x93) + x92) +
                                           grad_test[2] * (x10 * (x16 * x89 + x96) + x94) + x43 * (x15 * x89 + x91)) +
                          grad_trial[2] * (grad_test[0] * (x10 * (x100 + x12 * x97) + x85) +
                                           grad_test[1] * (x10 * (x102 + x15 * x97) + x101) + x46 * (x16 * x97 + x99)));
                bf[3] +=
                    dx * (grad_trial[0] * (grad_test[1] * (x10 * (x22 * x37 + x68) + x55) +
                                           grad_test[2] * (x10 * (x25 * x37 + x77) + x76) + x38 * (x19 * x37 + x54)) +
                          grad_trial[1] * (grad_test[0] * (x10 * (x19 * x47 + x59) + x56) +
                                           grad_test[2] * (x10 * (x25 * x47 + x78) + x69) + x43 * (x22 * x47 + x67)) +
                          grad_trial[2] * (grad_test[0] * (x10 * (x19 * x51 + x64) + x60) +
                                           grad_test[1] * (x10 * (x22 * x51 + x72) + x70) + x46 * (x25 * x51 + x75)));
                bf[4] += dx * (grad_trial[0] * (x38 * (x11 + x19 * x52 + 4 * x19 * x53) + x43 * (x105 + x22 * x52) +
                                                x46 * (x106 + x25 * x52)) +
                               grad_trial[1] * (x38 * (x105 + x19 * x65) + x43 * (x11 + x22 * x65 + 4 * x22 * x66) +
                                                x46 * (x109 + x25 * x65)) +
                               grad_trial[2] * (x38 * (x106 + x19 * x73) + x43 * (x109 + x22 * x73) +
                                                x46 * (x11 + x25 * x73 + 4 * x25 * x74)));
                bf[5] +=
                    dx *
                    (grad_trial[0] * (grad_test[1] * (x10 * (x114 + x22 * x79) + x112) +
                                      grad_test[2] * (x10 * (x119 + x25 * x79) + x115) + x38 * (x110 + x19 * x79)) +
                     grad_trial[1] * (grad_test[0] * (x10 * (x121 + x19 * x89) + x111) +
                                      grad_test[2] * (x10 * (x125 + x25 * x89) + x123) + x43 * (x120 + x22 * x89)) +
                     grad_trial[2] * (grad_test[0] * (x10 * (x128 + x19 * x97) + x127) +
                                      grad_test[1] * (x10 * (x129 + x22 * x97) + x122) + x46 * (x126 + x25 * x97)));
                bf[6] +=
                    dx * (grad_trial[0] * (grad_test[1] * (x10 * (x31 * x37 + x93) + x92) +
                                           grad_test[2] * (x10 * (x100 + x34 * x37) + x85) + x38 * (x28 * x37 + x81)) +
                          grad_trial[1] * (grad_test[0] * (x10 * (x28 * x47 + x84) + x82) +
                                           grad_test[2] * (x10 * (x102 + x34 * x47) + x101) + x43 * (x31 * x47 + x91)) +
                          grad_trial[2] * (grad_test[0] * (x10 * (x28 * x51 + x88) + x86) +
                                           grad_test[1] * (x10 * (x31 * x51 + x96) + x94) + x46 * (x34 * x51 + x99)));
                bf[7] +=
                    dx *
                    (grad_trial[0] * (grad_test[1] * (x10 * (x121 + x31 * x52) + x111) +
                                      grad_test[2] * (x10 * (x128 + x34 * x52) + x127) + x38 * (x110 + x28 * x52)) +
                     grad_trial[1] * (grad_test[0] * (x10 * (x114 + x28 * x65) + x112) +
                                      grad_test[2] * (x10 * (x129 + x34 * x65) + x122) + x43 * (x120 + x31 * x65)) +
                     grad_trial[2] * (grad_test[0] * (x10 * (x119 + x28 * x73) + x115) +
                                      grad_test[1] * (x10 * (x125 + x31 * x73) + x123) + x46 * (x126 + x34 * x73)));
                bf[8] += dx * (grad_trial[0] * (x38 * (x11 + x28 * x79 + 4 * x28 * x80) + x43 * (x132 + x31 * x79) +
                                                x46 * (x133 + x34 * x79)) +
                               grad_trial[1] * (x38 * (x132 + x28 * x89) + x43 * (x11 + x31 * x89 + 4 * x31 * x90) +
                                                x46 * (x134 + x34 * x89)) +
                               grad_trial[2] * (x38 * (x133 + x28 * x97) + x43 * (x134 + x31 * x97) +
                                                x46 * (x11 + x34 * x97 + 4 * x34 * x98)));
            }

            T mu{1.0};
            T lambda{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_NeoHookeanWang_3_IMPL_hpp
