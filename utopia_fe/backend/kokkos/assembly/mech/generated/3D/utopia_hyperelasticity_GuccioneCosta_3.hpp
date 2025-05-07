#ifndef UTOPIA_TPL_HYPERELASTICITY_GuccioneCosta_3_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_GuccioneCosta_3_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_GuccioneCosta.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of GuccioneCosta for dimension 3
         */
        template <typename T>
        class GuccioneCosta<T, 3> {
        public:
            static constexpr int Dim = 3;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "GuccioneCosta_3"; }

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    in.get("mu", mu);
                    in.get("b_f", b_f);
                    in.get("b_t", b_t);
                    in.get("b_fs", b_fs);
                    in.get("k", k);
                }

                T mu{2000};
                T b_f{8};
                T b_t{2};
                T b_fs{4};
                T k{10};
            };

            GuccioneCosta(const Params &params) {
                mu = params.mu;
                b_f = params.b_f;
                b_t = params.b_t;
                b_fs = params.b_fs;
                k = params.k;
            }

            UTOPIA_FUNCTION void hessian(const T *UTOPIA_RESTRICT f,
                                         const T *grad_test,
                                         const T *grad_trial,
                                         const T dx,
                                         T *UTOPIA_RESTRICT bf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[3] * f[8];
                T x1 = 2 * f[5] * f[6] - 2 * x0;
                T x2 = k * (f[5] * f[6] - x0);
                T x3 = pow(f[0], 2);
                T x4 = b_fs * x3;
                T x5 = pow(f[1], 2);
                T x6 = pow(f[7], 2);
                T x7 = x6 - 1;
                T x8 = pow(f[2], 2);
                T x9 = pow(f[4], 2);
                T x10 = x8 + x9;
                T x11 = pow(f[3], 2);
                T x12 = pow(f[6], 2);
                T x13 = (1.0 / 2.0) * x11 + (1.0 / 2.0) * x12 + (1.0 / 2.0) * x3 - 1.0 / 2.0;
                T x14 = f[0] * f[1];
                T x15 = f[3] * f[4];
                T x16 = f[6] * f[7];
                T x17 = (1.0 / 2.0) * x14 + (1.0 / 2.0) * x15 + (1.0 / 2.0) * x16;
                T x18 = f[0] * f[2];
                T x19 = f[3] * f[5];
                T x20 = f[6] * f[8];
                T x21 = (1.0 / 2.0) * x18 + (1.0 / 2.0) * x19 + (1.0 / 2.0) * x20;
                T x22 = f[1] * f[2];
                T x23 = f[4] * f[5];
                T x24 = f[7] * f[8];
                T x25 = (1.0 / 2.0) * x22 + (1.0 / 2.0) * x23 + (1.0 / 2.0) * x24;
                T x26 = (1.0 / 2.0) * x5 + (1.0 / 2.0) * x6 + (1.0 / 2.0) * x9 - 1.0 / 2.0;
                T x27 = pow(f[5], 2);
                T x28 = pow(f[8], 2);
                T x29 = (1.0 / 2.0) * x27 + (1.0 / 2.0) * x28 + (1.0 / 2.0) * x8 - 1.0 / 2.0;
                T x30 = (1.0 / 2.0) * mu *
                        exp(b_f * pow(x13, 2) + b_fs * (2 * pow(x17, 2) + 2 * pow(x21, 2)) +
                            b_t * (2 * pow(x25, 2) + pow(x26, 2) + pow(x29, 2)));
                T x31 = 2 * x17;
                T x32 = b_fs * x31;
                T x33 = 2 * f[2];
                T x34 = 2 * f[1];
                T x35 = b_t * (x25 * x33 + x26 * x34) + f[0] * x32;
                T x36 = f[4] * f[8];
                T x37 = f[5] * f[7];
                T x38 = 2 * x36 - 2 * x37;
                T x39 = b_fs * x30;
                T x40 = 2 * b_f;
                T x41 = x13 * x40;
                T x42 = b_fs * (f[1] * x31 + x21 * x33) + f[0] * x41;
                T x43 = x30 * x42;
                T x44 = x35 * x43;
                T x45 = f[3] * f[7];
                T x46 = f[4] * f[6];
                T x47 = 2 * x45 - 2 * x46;
                T x48 = b_t * x30;
                T x49 = 2 * x21;
                T x50 = b_fs * x49;
                T x51 = b_t * (x25 * x34 + x29 * x33) + f[0] * x50;
                T x52 = x30 * x35;
                T x53 = x48 * (2 * x22 + x23 + x24) + x51 * x52;
                T x54 = k * (x45 - x46);
                T x55 = x27 + x5;
                T x56 = x28 - 1;
                T x57 = x43 * x51;
                T x58 = k * (x36 - x37);
                T x59 = f[1] * f[8];
                T x60 = 2 * f[2] * f[7] - 2 * x59;
                T x61 = f[0] * x40;
                T x62 = f[1] * f[4];
                T x63 = f[2] * f[5];
                T x64 = b_fs * (f[4] * x31 + f[5] * x49) + f[3] * x41;
                T x65 = x30 * (b_fs * (x62 + x63) + f[3] * x61) + x43 * x64;
                T x66 = f[0] * f[8];
                T x67 = f[2] * f[6];
                T x68 = 2 * x66 - 2 * x67;
                T x69 = 2 * f[8];
                T x70 = k * (f[0] * f[4] * f[8] - f[0] * x37 + f[1] * f[5] * f[6] - f[1] * x0 + f[2] * f[3] * f[7] -
                             f[2] * x46 - 1);
                T x71 = x69 * x70;
                T x72 = 2 * f[5];
                T x73 = 2 * f[4];
                T x74 = b_t * (x25 * x72 + x26 * x73) + f[3] * x32;
                T x75 = f[1] * f[3];
                T x76 = x39 * x75 + x43 * x74 + x71;
                T x77 = f[0] * f[7];
                T x78 = 2 * f[1] * f[6] - 2 * x77;
                T x79 = 2 * f[7];
                T x80 = x70 * x79;
                T x81 = b_t * (x25 * x73 + x29 * x72) + f[3] * x50;
                T x82 = f[2] * f[3] * x39 + x43 * x81 - x80;
                T x83 = b_fs * f[0];
                T x84 = f[3] * x83;
                T x85 = x30 * (b_t * (2 * x62 + x63) + x84) + x52 * x74;
                T x86 = f[0] * f[4];
                T x87 = x39 * x86 + x52 * x64 - x71;
                T x88 = 2 * x70;
                T x89 = f[6] * x88;
                T x90 = f[2] * f[4];
                T x91 = x48 * x90 + x52 * x81 + x89;
                T x92 = x30 * x51;
                T x93 = x30 * (b_t * (x62 + 2 * x63) + x84) + x81 * x92;
                T x94 = f[0] * f[5];
                T x95 = x39 * x94 + x64 * x92 + x80;
                T x96 = f[1] * f[5];
                T x97 = x48 * x96 + x74 * x92 - x89;
                T x98 = -2 * x90 + 2 * x96;
                T x99 = f[1] * f[7];
                T x100 = f[2] * f[8];
                T x101 = b_fs * (f[7] * x31 + f[8] * x49) + f[6] * x41;
                T x102 = x101 * x43 + x30 * (b_fs * (x100 + x99) + f[6] * x61);
                T x103 = 2 * f[2] * f[3] - 2 * x94;
                T x104 = x70 * x72;
                T x105 = b_t * (x25 * x69 + x26 * x79) + f[6] * x32;
                T x106 = f[1] * f[6] * x39 - x104 + x105 * x43;
                T x107 = -2 * x75 + 2 * x86;
                T x108 = x70 * x73;
                T x109 = b_t * (x25 * x79 + x29 * x69) + f[6] * x50;
                T x110 = x108 + x109 * x43 + x39 * x67;
                T x111 = f[6] * x83;
                T x112 = x105 * x52 + x30 * (b_t * (x100 + 2 * x99) + x111);
                T x113 = x101 * x52 + x104 + x39 * x77;
                T x114 = f[3] * x88;
                T x115 = f[2] * f[7] * x48 + x109 * x52 - x114;
                T x116 = x109 * x92 + x30 * (b_t * (2 * x100 + x99) + x111);
                T x117 = x101 * x92 - x108 + x39 * x66;
                T x118 = x105 * x92 + x114 + x48 * x59;
                T x119 = k * (f[2] * f[7] - x59);
                T x120 = k * (x66 - x67);
                T x121 = k * (f[1] * f[6] - x77);
                T x122 = b_fs * x11;
                T x123 = x30 * x64;
                T x124 = x123 * x74;
                T x125 = x30 * x74;
                T x126 = x125 * x81 + x48 * (x22 + 2 * x23 + x24);
                T x127 = x123 * x81;
                T x128 = f[3] * f[6];
                T x129 = f[4] * f[7];
                T x130 = f[5] * f[8];
                T x131 = x101 * x123 + x30 * (b_fs * (x129 + x130) + x128 * x40);
                T x132 = x33 * x70;
                T x133 = x105 * x123 + x132 + x39 * x46;
                T x134 = x34 * x70;
                T x135 = f[5] * f[6] * x39 + x109 * x123 - x134;
                T x136 = b_fs * x128;
                T x137 = x105 * x125 + x30 * (b_t * (2 * x129 + x130) + x136);
                T x138 = x101 * x125 - x132 + x39 * x45;
                T x139 = f[0] * x88;
                T x140 = x109 * x125 + x139 + x37 * x48;
                T x141 = x30 * x81;
                T x142 = x109 * x141 + x30 * (b_t * (x129 + 2 * x130) + x136);
                T x143 = x0 * x39 + x101 * x141 + x134;
                T x144 = x105 * x141 - x139 + x36 * x48;
                T x145 = k * (-x90 + x96);
                T x146 = k * (f[2] * f[3] - x94);
                T x147 = k * (-x75 + x86);
                T x148 = b_fs * x12;
                T x149 = x101 * x30;
                T x150 = x105 * x149;
                T x151 = x105 * x109 * x30 + x48 * (x22 + x23 + 2 * x24);
                T x152 = x109 * x149;
                bf[0] +=
                    dx *
                    (grad_test[0] *
                         (grad_trial[0] * (x30 * pow(x42, 2) + x30 * (b_fs * (x5 + x8) + x3 * x40 + x41) + x38 * x58) +
                          grad_trial[1] * (x1 * x58 + x30 * (b_fs * x14 + x32) + x44) +
                          grad_trial[2] * (x30 * (b_fs * x18 + x50) + x47 * x58 + x57)) +
                     grad_test[1] *
                         (grad_trial[0] * (x2 * x38 + x39 * (2 * x14 + x15 + x16) + x44) +
                          grad_trial[1] * (x1 * x2 + x30 * pow(x35, 2) + x30 * (b_t * (x10 + 3 * x5 + x7) + x4)) +
                          grad_trial[2] * (x2 * x47 + x53)) +
                     grad_test[2] *
                         (grad_trial[0] * (x38 * x54 + x39 * (2 * x18 + x19 + x20) + x57) +
                          grad_trial[1] * (x1 * x54 + x53) +
                          grad_trial[2] * (x30 * pow(x51, 2) + x30 * (b_t * (x55 + x56 + 3 * x8) + x4) + x47 * x54)));
                bf[1] += dx * (grad_test[0] * (grad_trial[0] * (x58 * x60 + x65) + grad_trial[1] * (x58 * x68 + x76) +
                                               grad_trial[2] * (x58 * x78 + x82)) +
                               grad_test[1] * (grad_trial[0] * (x2 * x60 + x87) + grad_trial[1] * (x2 * x68 + x85) +
                                               grad_trial[2] * (x2 * x78 + x91)) +
                               grad_test[2] * (grad_trial[0] * (x54 * x60 + x95) + grad_trial[1] * (x54 * x68 + x97) +
                                               grad_trial[2] * (x54 * x78 + x93)));
                bf[2] +=
                    dx * (grad_test[0] * (grad_trial[0] * (x102 + x58 * x98) + grad_trial[1] * (x103 * x58 + x106) +
                                          grad_trial[2] * (x107 * x58 + x110)) +
                          grad_test[1] * (grad_trial[0] * (x113 + x2 * x98) + grad_trial[1] * (x103 * x2 + x112) +
                                          grad_trial[2] * (x107 * x2 + x115)) +
                          grad_test[2] * (grad_trial[0] * (x117 + x54 * x98) + grad_trial[1] * (x103 * x54 + x118) +
                                          grad_trial[2] * (x107 * x54 + x116)));
                bf[3] += dx * (grad_test[0] * (grad_trial[0] * (x119 * x38 + x65) + grad_trial[1] * (x1 * x119 + x87) +
                                               grad_trial[2] * (x119 * x47 + x95)) +
                               grad_test[1] * (grad_trial[0] * (x120 * x38 + x76) + grad_trial[1] * (x1 * x120 + x85) +
                                               grad_trial[2] * (x120 * x47 + x97)) +
                               grad_test[2] * (grad_trial[0] * (x121 * x38 + x82) + grad_trial[1] * (x1 * x121 + x91) +
                                               grad_trial[2] * (x121 * x47 + x93)));
                bf[4] += dx * (grad_test[0] * (grad_trial[0] * (x119 * x60 + x30 * pow(x64, 2) +
                                                                x30 * (b_fs * (x27 + x9) + x11 * x40 + x41)) +
                                               grad_trial[1] * (x119 * x68 + x124 + x30 * (b_fs * x15 + x32)) +
                                               grad_trial[2] * (x119 * x78 + x127 + x30 * (b_fs * x19 + x50))) +
                               grad_test[1] * (grad_trial[0] * (x120 * x60 + x124 + x39 * (x14 + 2 * x15 + x16)) +
                                               grad_trial[1] * (x120 * x68 + x30 * pow(x74, 2) +
                                                                x30 * (b_t * (x55 + x7 + 3 * x9) + x122)) +
                                               grad_trial[2] * (x120 * x78 + x126)) +
                               grad_test[2] * (grad_trial[0] * (x121 * x60 + x127 + x39 * (x18 + 2 * x19 + x20)) +
                                               grad_trial[1] * (x121 * x68 + x126) +
                                               grad_trial[2] * (x121 * x78 + x30 * pow(x81, 2) +
                                                                x30 * (b_t * (x10 + 3 * x27 + x56) + x122))));
                bf[5] +=
                    dx * (grad_test[0] * (grad_trial[0] * (x119 * x98 + x131) + grad_trial[1] * (x103 * x119 + x133) +
                                          grad_trial[2] * (x107 * x119 + x135)) +
                          grad_test[1] * (grad_trial[0] * (x120 * x98 + x138) + grad_trial[1] * (x103 * x120 + x137) +
                                          grad_trial[2] * (x107 * x120 + x140)) +
                          grad_test[2] * (grad_trial[0] * (x121 * x98 + x143) + grad_trial[1] * (x103 * x121 + x144) +
                                          grad_trial[2] * (x107 * x121 + x142)));
                bf[6] +=
                    dx * (grad_test[0] * (grad_trial[0] * (x102 + x145 * x38) + grad_trial[1] * (x1 * x145 + x113) +
                                          grad_trial[2] * (x117 + x145 * x47)) +
                          grad_test[1] * (grad_trial[0] * (x106 + x146 * x38) + grad_trial[1] * (x1 * x146 + x112) +
                                          grad_trial[2] * (x118 + x146 * x47)) +
                          grad_test[2] * (grad_trial[0] * (x110 + x147 * x38) + grad_trial[1] * (x1 * x147 + x115) +
                                          grad_trial[2] * (x116 + x147 * x47)));
                bf[7] +=
                    dx * (grad_test[0] * (grad_trial[0] * (x131 + x145 * x60) + grad_trial[1] * (x138 + x145 * x68) +
                                          grad_trial[2] * (x143 + x145 * x78)) +
                          grad_test[1] * (grad_trial[0] * (x133 + x146 * x60) + grad_trial[1] * (x137 + x146 * x68) +
                                          grad_trial[2] * (x144 + x146 * x78)) +
                          grad_test[2] * (grad_trial[0] * (x135 + x147 * x60) + grad_trial[1] * (x140 + x147 * x68) +
                                          grad_trial[2] * (x142 + x147 * x78)));
                bf[8] += dx * (grad_test[0] * (grad_trial[0] * (pow(x101, 2) * x30 + x145 * x98 +
                                                                x30 * (b_fs * (x28 + x6) + x12 * x40 + x41)) +
                                               grad_trial[1] * (x103 * x145 + x150 + x30 * (b_fs * x16 + x32)) +
                                               grad_trial[2] * (x107 * x145 + x152 + x30 * (b_fs * x20 + x50))) +
                               grad_test[1] * (grad_trial[0] * (x146 * x98 + x150 + x39 * (x14 + x15 + 2 * x16)) +
                                               grad_trial[1] * (x103 * x146 + pow(x105, 2) * x30 +
                                                                x30 * (b_t * (x5 + x56 + 3 * x6 + x9) + x148)) +
                                               grad_trial[2] * (x107 * x146 + x151)) +
                               grad_test[2] * (grad_trial[0] * (x147 * x98 + x152 + x39 * (x18 + x19 + 2 * x20)) +
                                               grad_trial[1] * (x103 * x147 + x151) +
                                               grad_trial[2] * (x107 * x147 + pow(x109, 2) * x30 +
                                                                x30 * (b_t * (x27 + 3 * x28 + x7 + x8) + x148))));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[5] * f[7];
                T x1 = f[3] * f[8];
                T x2 = f[4] * f[6];
                T x3 = k * (f[0] * f[4] * f[8] - f[0] * x0 + f[1] * f[5] * f[6] - f[1] * x1 + f[2] * f[3] * f[7] -
                            f[2] * x2 - 1);
                T x4 = 2 * f[0];
                T x5 = (1.0 / 2.0) * pow(f[0], 2) + (1.0 / 2.0) * pow(f[3], 2) + (1.0 / 2.0) * pow(f[6], 2) - 1.0 / 2.0;
                T x6 = b_f * x5;
                T x7 = (1.0 / 2.0) * f[0];
                T x8 = (1.0 / 2.0) * f[3];
                T x9 = (1.0 / 2.0) * f[6];
                T x10 = f[1] * x7 + f[4] * x8 + f[7] * x9;
                T x11 = 2 * f[1];
                T x12 = f[2] * x7 + f[5] * x8 + f[8] * x9;
                T x13 = 2 * f[2];
                T x14 = (1.0 / 2.0) * f[1] * f[2] + (1.0 / 2.0) * f[4] * f[5] + (1.0 / 2.0) * f[7] * f[8];
                T x15 =
                    (1.0 / 2.0) * pow(f[1], 2) + (1.0 / 2.0) * pow(f[4], 2) + (1.0 / 2.0) * pow(f[7], 2) - 1.0 / 2.0;
                T x16 =
                    (1.0 / 2.0) * pow(f[2], 2) + (1.0 / 2.0) * pow(f[5], 2) + (1.0 / 2.0) * pow(f[8], 2) - 1.0 / 2.0;
                T x17 = (1.0 / 2.0) * mu *
                        exp(b_f * pow(x5, 2) + b_fs * (2 * pow(x10, 2) + 2 * pow(x12, 2)) +
                            b_t * (2 * pow(x14, 2) + pow(x15, 2) + pow(x16, 2)));
                T x18 = b_fs * x4;
                T x19 = 2 * f[3];
                T x20 = 2 * f[4];
                T x21 = 2 * f[5];
                T x22 = b_fs * x19;
                T x23 = 2 * f[6];
                T x24 = 2 * f[7];
                T x25 = 2 * f[8];
                T x26 = b_fs * x23;
                lf[0] +=
                    dx * (grad_test[0] *
                              (x17 * (b_fs * (x10 * x11 + x12 * x13) + x4 * x6) + x3 * (2 * f[4] * f[8] - 2 * x0)) +
                          grad_test[1] *
                              (x17 * (b_t * (x11 * x15 + x13 * x14) + x10 * x18) + x3 * (2 * f[5] * f[6] - 2 * x1)) +
                          grad_test[2] *
                              (x17 * (b_t * (x11 * x14 + x13 * x16) + x12 * x18) + x3 * (2 * f[3] * f[7] - 2 * x2)));
                lf[1] +=
                    dx * (grad_test[0] * (x17 * (b_fs * (x10 * x20 + x12 * x21) + x19 * x6) +
                                          x3 * (2 * f[2] * f[7] - f[8] * x11)) +
                          grad_test[1] *
                              (x17 * (b_t * (x14 * x21 + x15 * x20) + x10 * x22) + x3 * (-f[6] * x13 + f[8] * x4)) +
                          grad_test[2] *
                              (x17 * (b_t * (x14 * x20 + x16 * x21) + x12 * x22) + x3 * (2 * f[1] * f[6] - f[7] * x4)));
                lf[2] +=
                    dx * (grad_test[0] *
                              (x17 * (b_fs * (x10 * x24 + x12 * x25) + x23 * x6) + x3 * (-f[4] * x13 + f[5] * x11)) +
                          grad_test[1] * (x17 * (b_t * (x14 * x25 + x15 * x24) + x10 * x26) +
                                          x3 * (-f[0] * x21 + 2 * f[2] * f[3])) +
                          grad_test[2] *
                              (x17 * (b_t * (x14 * x24 + x16 * x25) + x12 * x26) + x3 * (f[0] * x20 - f[3] * x11)));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = (1.0 / 2.0) * f[0];
                T x1 = (1.0 / 2.0) * f[3];
                T x2 = (1.0 / 2.0) * f[6];
                e += dx * (k * pow(f[0] * f[4] * f[8] - f[0] * f[5] * f[7] - f[1] * f[3] * f[8] + f[1] * f[5] * f[6] +
                                       f[2] * f[3] * f[7] - f[2] * f[4] * f[6] - 1,
                                   2) +
                           (1.0 / 2.0) * mu *
                               (exp(b_f * pow((1.0 / 2.0) * pow(f[0], 2) + (1.0 / 2.0) * pow(f[3], 2) +
                                                  (1.0 / 2.0) * pow(f[6], 2) - 1.0 / 2.0,
                                              2) +
                                    b_fs * (2 * pow(f[1] * x0 + f[4] * x1 + f[7] * x2, 2) +
                                            2 * pow(f[2] * x0 + f[5] * x1 + f[8] * x2, 2)) +
                                    b_t * (2 * pow((1.0 / 2.0) * f[1] * f[2] + (1.0 / 2.0) * f[4] * f[5] +
                                                       (1.0 / 2.0) * f[7] * f[8],
                                                   2) +
                                           pow((1.0 / 2.0) * pow(f[1], 2) + (1.0 / 2.0) * pow(f[4], 2) +
                                                   (1.0 / 2.0) * pow(f[7], 2) - 1.0 / 2.0,
                                               2) +
                                           pow((1.0 / 2.0) * pow(f[2], 2) + (1.0 / 2.0) * pow(f[5], 2) +
                                                   (1.0 / 2.0) * pow(f[8], 2) - 1.0 / 2.0,
                                               2))) -
                                1));
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
                T x0 = f[5] * f[7];
                T x1 = f[3] * f[8];
                T x2 = f[4] * f[6];
                T x3 = f[0] * f[4] * f[8] - f[0] * x0 + f[1] * f[5] * f[6] - f[1] * x1 + f[2] * f[3] * f[7] -
                       f[2] * x2 - 1;
                T x4 = pow(f[0], 2);
                T x5 = pow(f[3], 2);
                T x6 = pow(f[6], 2);
                T x7 = (1.0 / 2.0) * x4 + (1.0 / 2.0) * x5 + (1.0 / 2.0) * x6 - 1.0 / 2.0;
                T x8 = f[0] * f[1];
                T x9 = f[3] * f[4];
                T x10 = f[6] * f[7];
                T x11 = (1.0 / 2.0) * x10 + (1.0 / 2.0) * x8 + (1.0 / 2.0) * x9;
                T x12 = f[0] * f[2];
                T x13 = f[3] * f[5];
                T x14 = f[6] * f[8];
                T x15 = (1.0 / 2.0) * x12 + (1.0 / 2.0) * x13 + (1.0 / 2.0) * x14;
                T x16 = f[1] * f[2];
                T x17 = f[4] * f[5];
                T x18 = f[7] * f[8];
                T x19 = (1.0 / 2.0) * x16 + (1.0 / 2.0) * x17 + (1.0 / 2.0) * x18;
                T x20 = pow(f[1], 2);
                T x21 = pow(f[4], 2);
                T x22 = pow(f[7], 2);
                T x23 = (1.0 / 2.0) * x20 + (1.0 / 2.0) * x21 + (1.0 / 2.0) * x22 - 1.0 / 2.0;
                T x24 = pow(f[2], 2);
                T x25 = pow(f[5], 2);
                T x26 = pow(f[8], 2);
                T x27 = (1.0 / 2.0) * x24 + (1.0 / 2.0) * x25 + (1.0 / 2.0) * x26 - 1.0 / 2.0;
                T x28 = exp(b_f * pow(x7, 2) + b_fs * (2 * pow(x11, 2) + 2 * pow(x15, 2)) +
                            b_t * (2 * pow(x19, 2) + pow(x23, 2) + pow(x27, 2)));
                T x29 = (1.0 / 2.0) * mu;
                T x30 = f[4] * f[8];
                T x31 = -2 * x0 + 2 * x30;
                T x32 = k * x3;
                T x33 = 2 * b_f;
                T x34 = x33 * x7;
                T x35 = 2 * f[1];
                T x36 = 2 * f[2];
                T x37 = b_fs * (x11 * x35 + x15 * x36) + f[0] * x34;
                T x38 = x28 * x29;
                T x39 = x37 * x38;
                T x40 = 2 * f[5] * f[6] - 2 * x1;
                T x41 = 2 * b_fs;
                T x42 = x11 * x41;
                T x43 = b_t * (x19 * x36 + x23 * x35) + f[0] * x42;
                T x44 = x38 * x43;
                T x45 = f[3] * f[7];
                T x46 = -2 * x2 + 2 * x45;
                T x47 = x15 * x41;
                T x48 = b_t * (x19 * x35 + x27 * x36) + f[0] * x47;
                T x49 = x38 * x48;
                T x50 = f[1] * f[8];
                T x51 = 2 * f[2] * f[7] - 2 * x50;
                T x52 = 2 * f[4];
                T x53 = 2 * f[5];
                T x54 = b_fs * (x11 * x52 + x15 * x53) + f[3] * x34;
                T x55 = x38 * x54;
                T x56 = f[0] * f[8];
                T x57 = f[2] * f[6];
                T x58 = 2 * x56 - 2 * x57;
                T x59 = b_t * (x19 * x53 + x23 * x52) + f[3] * x42;
                T x60 = x38 * x59;
                T x61 = f[0] * f[7];
                T x62 = 2 * f[1] * f[6] - 2 * x61;
                T x63 = b_t * (x19 * x52 + x27 * x53) + f[3] * x47;
                T x64 = x38 * x63;
                T x65 = f[1] * f[5];
                T x66 = f[2] * f[4];
                T x67 = 2 * x65 - 2 * x66;
                T x68 = 2 * f[7];
                T x69 = 2 * f[8];
                T x70 = b_fs * (x11 * x68 + x15 * x69) + f[6] * x34;
                T x71 = x38 * x70;
                T x72 = f[0] * f[5];
                T x73 = 2 * f[2] * f[3] - 2 * x72;
                T x74 = b_t * (x19 * x69 + x23 * x68) + f[6] * x42;
                T x75 = x38 * x74;
                T x76 = f[0] * f[4];
                T x77 = f[1] * f[3];
                T x78 = 2 * x76 - 2 * x77;
                T x79 = b_t * (x19 * x68 + x27 * x69) + f[6] * x47;
                T x80 = k * (f[5] * f[6] - x1);
                T x81 = b_fs * x4;
                T x82 = x22 - 1;
                T x83 = x21 + x24;
                T x84 = b_fs * x38;
                T x85 = x39 * x43;
                T x86 = b_t * x38;
                T x87 = x44 * x48 + x86 * (2 * x16 + x17 + x18);
                T x88 = k * (-x2 + x45);
                T x89 = x20 + x25;
                T x90 = x26 - 1;
                T x91 = x39 * x48;
                T x92 = k * (-x0 + x30);
                T x93 = f[0] * x33;
                T x94 = f[1] * f[4];
                T x95 = f[2] * f[5];
                T x96 = x38 * (b_fs * (x94 + x95) + f[3] * x93) + x39 * x54;
                T x97 = x32 * x69;
                T x98 = x39 * x59 + x77 * x84 + x97;
                T x99 = x32 * x68;
                T x100 = f[2] * f[3] * x84 + x39 * x63 - x99;
                T x101 = b_fs * f[0];
                T x102 = f[3] * x101;
                T x103 = x38 * (b_t * (2 * x94 + x95) + x102) + x44 * x59;
                T x104 = x44 * x54 + x76 * x84 - x97;
                T x105 = 2 * x32;
                T x106 = f[6] * x105;
                T x107 = x106 + x44 * x63 + x66 * x86;
                T x108 = x38 * (b_t * (x94 + 2 * x95) + x102) + x49 * x63;
                T x109 = x49 * x54 + x72 * x84 + x99;
                T x110 = -x106 + x49 * x59 + x65 * x86;
                T x111 = f[1] * f[7];
                T x112 = f[2] * f[8];
                T x113 = x38 * (b_fs * (x111 + x112) + f[6] * x93) + x39 * x70;
                T x114 = x32 * x53;
                T x115 = f[1] * f[6] * x84 - x114 + x39 * x74;
                T x116 = x32 * x52;
                T x117 = x116 + x39 * x79 + x57 * x84;
                T x118 = f[6] * x101;
                T x119 = x38 * (b_t * (2 * x111 + x112) + x118) + x44 * x74;
                T x120 = x114 + x44 * x70 + x61 * x84;
                T x121 = f[3] * x105;
                T x122 = f[2] * f[7] * x86 - x121 + x44 * x79;
                T x123 = x38 * (b_t * (x111 + 2 * x112) + x118) + x49 * x79;
                T x124 = -x116 + x49 * x70 + x56 * x84;
                T x125 = x121 + x49 * x74 + x50 * x86;
                T x126 = k * (f[2] * f[7] - x50);
                T x127 = k * (x56 - x57);
                T x128 = k * (f[1] * f[6] - x61);
                T x129 = b_fs * x5;
                T x130 = x55 * x59;
                T x131 = x60 * x63 + x86 * (x16 + 2 * x17 + x18);
                T x132 = x55 * x63;
                T x133 = f[3] * f[6];
                T x134 = f[4] * f[7];
                T x135 = f[5] * f[8];
                T x136 = x38 * (b_fs * (x134 + x135) + x133 * x33) + x55 * x70;
                T x137 = x32 * x36;
                T x138 = x137 + x2 * x84 + x55 * x74;
                T x139 = x32 * x35;
                T x140 = f[5] * f[6] * x84 - x139 + x55 * x79;
                T x141 = b_fs * x133;
                T x142 = x38 * (b_t * (2 * x134 + x135) + x141) + x60 * x74;
                T x143 = -x137 + x45 * x84 + x60 * x70;
                T x144 = f[0] * x105;
                T x145 = x0 * x86 + x144 + x60 * x79;
                T x146 = x38 * (b_t * (x134 + 2 * x135) + x141) + x64 * x79;
                T x147 = x1 * x84 + x139 + x64 * x70;
                T x148 = -x144 + x30 * x86 + x64 * x74;
                T x149 = k * (x65 - x66);
                T x150 = k * (f[2] * f[3] - x72);
                T x151 = k * (x76 - x77);
                T x152 = b_fs * x6;
                T x153 = x71 * x74;
                T x154 = x75 * x79 + x86 * (x16 + x17 + 2 * x18);
                T x155 = x71 * x79;
                e += dx * (k * pow(x3, 2) + x29 * (x28 - 1));
                lf[0] += dx * (grad_test[0] * (x31 * x32 + x39) + grad_test[1] * (x32 * x40 + x44) +
                               grad_test[2] * (x32 * x46 + x49));
                lf[1] += dx * (grad_test[0] * (x32 * x51 + x55) + grad_test[1] * (x32 * x58 + x60) +
                               grad_test[2] * (x32 * x62 + x64));
                lf[2] += dx * (grad_test[0] * (x32 * x67 + x71) + grad_test[1] * (x32 * x73 + x75) +
                               grad_test[2] * (x32 * x78 + x38 * x79));
                bf[0] +=
                    dx *
                    (grad_test[0] * (grad_trial[0] *
                                         (x31 * x92 + pow(x37, 2) * x38 + x38 * (b_fs * (x20 + x24) + x33 * x4 + x34)) +
                                     grad_trial[1] * (x38 * (b_fs * x8 + x42) + x40 * x92 + x85) +
                                     grad_trial[2] * (x38 * (b_fs * x12 + x47) + x46 * x92 + x91)) +
                     grad_test[1] *
                         (grad_trial[0] * (x31 * x80 + x84 * (x10 + 2 * x8 + x9) + x85) +
                          grad_trial[1] * (x38 * pow(x43, 2) + x38 * (b_t * (3 * x20 + x82 + x83) + x81) + x40 * x80) +
                          grad_trial[2] * (x46 * x80 + x87)) +
                     grad_test[2] *
                         (grad_trial[0] * (x31 * x88 + x84 * (2 * x12 + x13 + x14) + x91) +
                          grad_trial[1] * (x40 * x88 + x87) +
                          grad_trial[2] * (x38 * pow(x48, 2) + x38 * (b_t * (3 * x24 + x89 + x90) + x81) + x46 * x88)));
                bf[1] += dx * (grad_test[0] * (grad_trial[0] * (x51 * x92 + x96) + grad_trial[1] * (x58 * x92 + x98) +
                                               grad_trial[2] * (x100 + x62 * x92)) +
                               grad_test[1] * (grad_trial[0] * (x104 + x51 * x80) + grad_trial[1] * (x103 + x58 * x80) +
                                               grad_trial[2] * (x107 + x62 * x80)) +
                               grad_test[2] * (grad_trial[0] * (x109 + x51 * x88) + grad_trial[1] * (x110 + x58 * x88) +
                                               grad_trial[2] * (x108 + x62 * x88)));
                bf[2] += dx * (grad_test[0] * (grad_trial[0] * (x113 + x67 * x92) + grad_trial[1] * (x115 + x73 * x92) +
                                               grad_trial[2] * (x117 + x78 * x92)) +
                               grad_test[1] * (grad_trial[0] * (x120 + x67 * x80) + grad_trial[1] * (x119 + x73 * x80) +
                                               grad_trial[2] * (x122 + x78 * x80)) +
                               grad_test[2] * (grad_trial[0] * (x124 + x67 * x88) + grad_trial[1] * (x125 + x73 * x88) +
                                               grad_trial[2] * (x123 + x78 * x88)));
                bf[3] +=
                    dx * (grad_test[0] * (grad_trial[0] * (x126 * x31 + x96) + grad_trial[1] * (x104 + x126 * x40) +
                                          grad_trial[2] * (x109 + x126 * x46)) +
                          grad_test[1] * (grad_trial[0] * (x127 * x31 + x98) + grad_trial[1] * (x103 + x127 * x40) +
                                          grad_trial[2] * (x110 + x127 * x46)) +
                          grad_test[2] * (grad_trial[0] * (x100 + x128 * x31) + grad_trial[1] * (x107 + x128 * x40) +
                                          grad_trial[2] * (x108 + x128 * x46)));
                bf[4] += dx * (grad_test[0] * (grad_trial[0] * (x126 * x51 + x38 * pow(x54, 2) +
                                                                x38 * (b_fs * (x21 + x25) + x33 * x5 + x34)) +
                                               grad_trial[1] * (x126 * x58 + x130 + x38 * (b_fs * x9 + x42)) +
                                               grad_trial[2] * (x126 * x62 + x132 + x38 * (b_fs * x13 + x47))) +
                               grad_test[1] * (grad_trial[0] * (x127 * x51 + x130 + x84 * (x10 + x8 + 2 * x9)) +
                                               grad_trial[1] * (x127 * x58 + x38 * pow(x59, 2) +
                                                                x38 * (b_t * (3 * x21 + x82 + x89) + x129)) +
                                               grad_trial[2] * (x127 * x62 + x131)) +
                               grad_test[2] * (grad_trial[0] * (x128 * x51 + x132 + x84 * (x12 + 2 * x13 + x14)) +
                                               grad_trial[1] * (x128 * x58 + x131) +
                                               grad_trial[2] * (x128 * x62 + x38 * pow(x63, 2) +
                                                                x38 * (b_t * (3 * x25 + x83 + x90) + x129))));
                bf[5] +=
                    dx * (grad_test[0] * (grad_trial[0] * (x126 * x67 + x136) + grad_trial[1] * (x126 * x73 + x138) +
                                          grad_trial[2] * (x126 * x78 + x140)) +
                          grad_test[1] * (grad_trial[0] * (x127 * x67 + x143) + grad_trial[1] * (x127 * x73 + x142) +
                                          grad_trial[2] * (x127 * x78 + x145)) +
                          grad_test[2] * (grad_trial[0] * (x128 * x67 + x147) + grad_trial[1] * (x128 * x73 + x148) +
                                          grad_trial[2] * (x128 * x78 + x146)));
                bf[6] +=
                    dx * (grad_test[0] * (grad_trial[0] * (x113 + x149 * x31) + grad_trial[1] * (x120 + x149 * x40) +
                                          grad_trial[2] * (x124 + x149 * x46)) +
                          grad_test[1] * (grad_trial[0] * (x115 + x150 * x31) + grad_trial[1] * (x119 + x150 * x40) +
                                          grad_trial[2] * (x125 + x150 * x46)) +
                          grad_test[2] * (grad_trial[0] * (x117 + x151 * x31) + grad_trial[1] * (x122 + x151 * x40) +
                                          grad_trial[2] * (x123 + x151 * x46)));
                bf[7] +=
                    dx * (grad_test[0] * (grad_trial[0] * (x136 + x149 * x51) + grad_trial[1] * (x143 + x149 * x58) +
                                          grad_trial[2] * (x147 + x149 * x62)) +
                          grad_test[1] * (grad_trial[0] * (x138 + x150 * x51) + grad_trial[1] * (x142 + x150 * x58) +
                                          grad_trial[2] * (x148 + x150 * x62)) +
                          grad_test[2] * (grad_trial[0] * (x140 + x151 * x51) + grad_trial[1] * (x145 + x151 * x58) +
                                          grad_trial[2] * (x146 + x151 * x62)));
                bf[8] += dx * (grad_test[0] * (grad_trial[0] * (x149 * x67 + x38 * pow(x70, 2) +
                                                                x38 * (b_fs * (x22 + x26) + x33 * x6 + x34)) +
                                               grad_trial[1] * (x149 * x73 + x153 + x38 * (b_fs * x10 + x42)) +
                                               grad_trial[2] * (x149 * x78 + x155 + x38 * (b_fs * x14 + x47))) +
                               grad_test[1] * (grad_trial[0] * (x150 * x67 + x153 + x84 * (2 * x10 + x8 + x9)) +
                                               grad_trial[1] * (x150 * x73 + x38 * pow(x74, 2) +
                                                                x38 * (b_t * (x20 + x21 + 3 * x22 + x90) + x152)) +
                                               grad_trial[2] * (x150 * x78 + x154)) +
                               grad_test[2] * (grad_trial[0] * (x151 * x67 + x155 + x84 * (x12 + x13 + 2 * x14)) +
                                               grad_trial[1] * (x151 * x73 + x154) +
                                               grad_trial[2] * (x151 * x78 + x38 * pow(x79, 2) +
                                                                x38 * (b_t * (x24 + x25 + 3 * x26 + x82) + x152))));
            }

            UTOPIA_FUNCTION void apply(const T *UTOPIA_RESTRICT f,
                                       const T *grad_test,
                                       const T *disp_grad,
                                       const T dx,
                                       T *UTOPIA_RESTRICT res) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[3] * f[8];
                T x1 = 2 * f[5] * f[6] - 2 * x0;
                T x2 = k * (f[5] * f[6] - x0);
                T x3 = pow(f[0], 2);
                T x4 = b_fs * x3;
                T x5 = pow(f[1], 2);
                T x6 = pow(f[7], 2);
                T x7 = x6 - 1;
                T x8 = pow(f[2], 2);
                T x9 = pow(f[4], 2);
                T x10 = x8 + x9;
                T x11 = pow(f[3], 2);
                T x12 = pow(f[6], 2);
                T x13 = (1.0 / 2.0) * x11 + (1.0 / 2.0) * x12 + (1.0 / 2.0) * x3 - 1.0 / 2.0;
                T x14 = f[0] * f[1];
                T x15 = f[3] * f[4];
                T x16 = f[6] * f[7];
                T x17 = (1.0 / 2.0) * x14 + (1.0 / 2.0) * x15 + (1.0 / 2.0) * x16;
                T x18 = f[0] * f[2];
                T x19 = f[3] * f[5];
                T x20 = f[6] * f[8];
                T x21 = (1.0 / 2.0) * x18 + (1.0 / 2.0) * x19 + (1.0 / 2.0) * x20;
                T x22 = f[1] * f[2];
                T x23 = f[4] * f[5];
                T x24 = f[7] * f[8];
                T x25 = (1.0 / 2.0) * x22 + (1.0 / 2.0) * x23 + (1.0 / 2.0) * x24;
                T x26 = (1.0 / 2.0) * x5 + (1.0 / 2.0) * x6 + (1.0 / 2.0) * x9 - 1.0 / 2.0;
                T x27 = pow(f[5], 2);
                T x28 = pow(f[8], 2);
                T x29 = (1.0 / 2.0) * x27 + (1.0 / 2.0) * x28 + (1.0 / 2.0) * x8 - 1.0 / 2.0;
                T x30 = (1.0 / 2.0) * mu *
                        exp(b_f * pow(x13, 2) + b_fs * (2 * pow(x17, 2) + 2 * pow(x21, 2)) +
                            b_t * (2 * pow(x25, 2) + pow(x26, 2) + pow(x29, 2)));
                T x31 = 2 * x17;
                T x32 = b_fs * x31;
                T x33 = 2 * f[2];
                T x34 = 2 * f[1];
                T x35 = b_t * (x25 * x33 + x26 * x34) + f[0] * x32;
                T x36 = f[4] * f[8];
                T x37 = f[5] * f[7];
                T x38 = 2 * x36 - 2 * x37;
                T x39 = b_fs * x30;
                T x40 = 2 * b_f;
                T x41 = x13 * x40;
                T x42 = b_fs * (f[1] * x31 + x21 * x33) + f[0] * x41;
                T x43 = x30 * x42;
                T x44 = x35 * x43;
                T x45 = f[3] * f[7];
                T x46 = f[4] * f[6];
                T x47 = 2 * x45 - 2 * x46;
                T x48 = b_t * x30;
                T x49 = 2 * x21;
                T x50 = b_fs * x49;
                T x51 = b_t * (x25 * x34 + x29 * x33) + f[0] * x50;
                T x52 = x30 * x35;
                T x53 = x48 * (2 * x22 + x23 + x24) + x51 * x52;
                T x54 = f[0] * f[8];
                T x55 = f[2] * f[6];
                T x56 = 2 * x54 - 2 * x55;
                T x57 = b_fs * f[0];
                T x58 = f[3] * x57;
                T x59 = f[2] * f[5];
                T x60 = f[1] * f[4];
                T x61 = 2 * f[5];
                T x62 = 2 * f[4];
                T x63 = b_t * (x25 * x61 + x26 * x62) + f[3] * x32;
                T x64 = x30 * (b_t * (x59 + 2 * x60) + x58) + x52 * x63;
                T x65 = f[0] * f[5];
                T x66 = 2 * f[2] * f[3] - 2 * x65;
                T x67 = f[6] * x57;
                T x68 = f[2] * f[8];
                T x69 = f[1] * f[7];
                T x70 = 2 * f[8];
                T x71 = 2 * f[7];
                T x72 = b_t * (x25 * x70 + x26 * x71) + f[6] * x32;
                T x73 = x30 * (b_t * (x68 + 2 * x69) + x67) + x52 * x72;
                T x74 = f[1] * f[8];
                T x75 = 2 * f[2] * f[7] - 2 * x74;
                T x76 = k * (f[0] * f[4] * f[8] - f[0] * x37 + f[1] * f[5] * f[6] - f[1] * x0 + f[2] * f[3] * f[7] -
                             f[2] * x46 - 1);
                T x77 = x70 * x76;
                T x78 = b_fs * (f[4] * x31 + f[5] * x49) + f[3] * x41;
                T x79 = f[0] * f[4];
                T x80 = x39 * x79 + x52 * x78 - x77;
                T x81 = f[0] * f[7];
                T x82 = 2 * f[1] * f[6] - 2 * x81;
                T x83 = 2 * x76;
                T x84 = f[6] * x83;
                T x85 = b_t * (x25 * x62 + x29 * x61) + f[3] * x50;
                T x86 = f[2] * f[4];
                T x87 = x48 * x86 + x52 * x85 + x84;
                T x88 = f[1] * f[5];
                T x89 = -2 * x86 + 2 * x88;
                T x90 = x61 * x76;
                T x91 = b_fs * (f[7] * x31 + f[8] * x49) + f[6] * x41;
                T x92 = x39 * x81 + x52 * x91 + x90;
                T x93 = f[1] * f[3];
                T x94 = 2 * x79 - 2 * x93;
                T x95 = f[3] * x83;
                T x96 = b_t * (x25 * x71 + x29 * x70) + f[6] * x50;
                T x97 = f[2] * f[7] * x48 + x52 * x96 - x95;
                T x98 = k * (x45 - x46);
                T x99 = x27 + x5;
                T x100 = x28 - 1;
                T x101 = x43 * x51;
                T x102 = x30 * x51;
                T x103 = x102 * x85 + x30 * (b_t * (2 * x59 + x60) + x58);
                T x104 = x102 * x96 + x30 * (b_t * (2 * x68 + x69) + x67);
                T x105 = x71 * x76;
                T x106 = x102 * x78 + x105 + x39 * x65;
                T x107 = x102 * x63 + x48 * x88 - x84;
                T x108 = x62 * x76;
                T x109 = x102 * x91 - x108 + x39 * x54;
                T x110 = x102 * x72 + x48 * x74 + x95;
                T x111 = k * (x36 - x37);
                T x112 = f[0] * x40;
                T x113 = x30 * (b_fs * (x59 + x60) + f[3] * x112) + x43 * x78;
                T x114 = x30 * (b_fs * (x68 + x69) + f[6] * x112) + x43 * x91;
                T x115 = x39 * x93 + x43 * x63 + x77;
                T x116 = f[2] * f[3] * x39 - x105 + x43 * x85;
                T x117 = f[1] * f[6] * x39 + x43 * x72 - x90;
                T x118 = x108 + x39 * x55 + x43 * x96;
                T x119 = k * (x54 - x55);
                T x120 = b_fs * x11;
                T x121 = x30 * x78;
                T x122 = x121 * x63;
                T x123 = x30 * x63;
                T x124 = x123 * x85 + x48 * (x22 + 2 * x23 + x24);
                T x125 = f[3] * f[6];
                T x126 = b_fs * x125;
                T x127 = f[5] * f[8];
                T x128 = f[4] * f[7];
                T x129 = x123 * x72 + x30 * (b_t * (x127 + 2 * x128) + x126);
                T x130 = x33 * x76;
                T x131 = x123 * x91 - x130 + x39 * x45;
                T x132 = f[0] * x83;
                T x133 = x123 * x96 + x132 + x37 * x48;
                T x134 = k * (f[1] * f[6] - x81);
                T x135 = x121 * x85;
                T x136 = x30 * x85;
                T x137 = x136 * x96 + x30 * (b_t * (2 * x127 + x128) + x126);
                T x138 = x34 * x76;
                T x139 = x0 * x39 + x136 * x91 + x138;
                T x140 = -x132 + x136 * x72 + x36 * x48;
                T x141 = k * (f[2] * f[7] - x74);
                T x142 = x121 * x91 + x30 * (b_fs * (x127 + x128) + x125 * x40);
                T x143 = x121 * x72 + x130 + x39 * x46;
                T x144 = f[5] * f[6] * x39 + x121 * x96 - x138;
                T x145 = k * (f[2] * f[3] - x65);
                T x146 = b_fs * x12;
                T x147 = x30 * x91;
                T x148 = x147 * x72;
                T x149 = x30 * x72 * x96 + x48 * (x22 + x23 + 2 * x24);
                T x150 = k * (x79 - x93);
                T x151 = x147 * x96;
                T x152 = k * (-x86 + x88);
                res[0] +=
                    dx *
                    (grad_test[0] *
                         (disp_grad[0] * (x111 * x38 + x30 * pow(x42, 2) + x30 * (b_fs * (x5 + x8) + x3 * x40 + x41)) +
                          disp_grad[1] * (x1 * x111 + x30 * (b_fs * x14 + x32) + x44) +
                          disp_grad[2] * (x101 + x111 * x47 + x30 * (b_fs * x18 + x50)) +
                          disp_grad[3] * (x111 * x75 + x113) + disp_grad[4] * (x111 * x56 + x115) +
                          disp_grad[5] * (x111 * x82 + x116) + disp_grad[6] * (x111 * x89 + x114) +
                          disp_grad[7] * (x111 * x66 + x117) + disp_grad[8] * (x111 * x94 + x118)) +
                     grad_test[1] *
                         (disp_grad[0] * (x2 * x38 + x39 * (2 * x14 + x15 + x16) + x44) +
                          disp_grad[1] * (x1 * x2 + x30 * pow(x35, 2) + x30 * (b_t * (x10 + 3 * x5 + x7) + x4)) +
                          disp_grad[2] * (x2 * x47 + x53) + disp_grad[3] * (x2 * x75 + x80) +
                          disp_grad[4] * (x2 * x56 + x64) + disp_grad[5] * (x2 * x82 + x87) +
                          disp_grad[6] * (x2 * x89 + x92) + disp_grad[7] * (x2 * x66 + x73) +
                          disp_grad[8] * (x2 * x94 + x97)) +
                     grad_test[2] *
                         (disp_grad[0] * (x101 + x38 * x98 + x39 * (2 * x18 + x19 + x20)) +
                          disp_grad[1] * (x1 * x98 + x53) +
                          disp_grad[2] * (x30 * pow(x51, 2) + x30 * (b_t * (x100 + 3 * x8 + x99) + x4) + x47 * x98) +
                          disp_grad[3] * (x106 + x75 * x98) + disp_grad[4] * (x107 + x56 * x98) +
                          disp_grad[5] * (x103 + x82 * x98) + disp_grad[6] * (x109 + x89 * x98) +
                          disp_grad[7] * (x110 + x66 * x98) + disp_grad[8] * (x104 + x94 * x98)));
                res[1] +=
                    dx * (grad_test[0] * (disp_grad[0] * (x113 + x141 * x38) + disp_grad[1] * (x1 * x141 + x80) +
                                          disp_grad[2] * (x106 + x141 * x47) +
                                          disp_grad[3] * (x141 * x75 + x30 * pow(x78, 2) +
                                                          x30 * (b_fs * (x27 + x9) + x11 * x40 + x41)) +
                                          disp_grad[4] * (x122 + x141 * x56 + x30 * (b_fs * x15 + x32)) +
                                          disp_grad[5] * (x135 + x141 * x82 + x30 * (b_fs * x19 + x50)) +
                                          disp_grad[6] * (x141 * x89 + x142) + disp_grad[7] * (x141 * x66 + x143) +
                                          disp_grad[8] * (x141 * x94 + x144)) +
                          grad_test[1] * (disp_grad[0] * (x115 + x119 * x38) + disp_grad[1] * (x1 * x119 + x64) +
                                          disp_grad[2] * (x107 + x119 * x47) +
                                          disp_grad[3] * (x119 * x75 + x122 + x39 * (x14 + 2 * x15 + x16)) +
                                          disp_grad[4] * (x119 * x56 + x30 * pow(x63, 2) +
                                                          x30 * (b_t * (x7 + 3 * x9 + x99) + x120)) +
                                          disp_grad[5] * (x119 * x82 + x124) + disp_grad[6] * (x119 * x89 + x131) +
                                          disp_grad[7] * (x119 * x66 + x129) + disp_grad[8] * (x119 * x94 + x133)) +
                          grad_test[2] * (disp_grad[0] * (x116 + x134 * x38) + disp_grad[1] * (x1 * x134 + x87) +
                                          disp_grad[2] * (x103 + x134 * x47) +
                                          disp_grad[3] * (x134 * x75 + x135 + x39 * (x18 + 2 * x19 + x20)) +
                                          disp_grad[4] * (x124 + x134 * x56) +
                                          disp_grad[5] * (x134 * x82 + x30 * pow(x85, 2) +
                                                          x30 * (b_t * (x10 + x100 + 3 * x27) + x120)) +
                                          disp_grad[6] * (x134 * x89 + x139) + disp_grad[7] * (x134 * x66 + x140) +
                                          disp_grad[8] * (x134 * x94 + x137)));
                res[2] +=
                    dx * (grad_test[0] * (disp_grad[0] * (x114 + x152 * x38) + disp_grad[1] * (x1 * x152 + x92) +
                                          disp_grad[2] * (x109 + x152 * x47) + disp_grad[3] * (x142 + x152 * x75) +
                                          disp_grad[4] * (x131 + x152 * x56) + disp_grad[5] * (x139 + x152 * x82) +
                                          disp_grad[6] * (x152 * x89 + x30 * pow(x91, 2) +
                                                          x30 * (b_fs * (x28 + x6) + x12 * x40 + x41)) +
                                          disp_grad[7] * (x148 + x152 * x66 + x30 * (b_fs * x16 + x32)) +
                                          disp_grad[8] * (x151 + x152 * x94 + x30 * (b_fs * x20 + x50))) +
                          grad_test[1] * (disp_grad[0] * (x117 + x145 * x38) + disp_grad[1] * (x1 * x145 + x73) +
                                          disp_grad[2] * (x110 + x145 * x47) + disp_grad[3] * (x143 + x145 * x75) +
                                          disp_grad[4] * (x129 + x145 * x56) + disp_grad[5] * (x140 + x145 * x82) +
                                          disp_grad[6] * (x145 * x89 + x148 + x39 * (x14 + x15 + 2 * x16)) +
                                          disp_grad[7] * (x145 * x66 + x30 * pow(x72, 2) +
                                                          x30 * (b_t * (x100 + x5 + 3 * x6 + x9) + x146)) +
                                          disp_grad[8] * (x145 * x94 + x149)) +
                          grad_test[2] * (disp_grad[0] * (x118 + x150 * x38) + disp_grad[1] * (x1 * x150 + x97) +
                                          disp_grad[2] * (x104 + x150 * x47) + disp_grad[3] * (x144 + x150 * x75) +
                                          disp_grad[4] * (x133 + x150 * x56) + disp_grad[5] * (x137 + x150 * x82) +
                                          disp_grad[6] * (x150 * x89 + x151 + x39 * (x18 + x19 + 2 * x20)) +
                                          disp_grad[7] * (x149 + x150 * x66) +
                                          disp_grad[8] * (x150 * x94 + x30 * pow(x96, 2) +
                                                          x30 * (b_t * (x27 + 3 * x28 + x7 + x8) + x146))));
            }

            T mu{2000};
            T b_f{8};
            T b_t{2};
            T b_fs{4};
            T k{10};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_GuccioneCosta_3_IMPL_hpp
