#ifndef UTOPIA_TPL_HYPERELASTICITY_NeoHookeanSiguenza_3_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeoHookeanSiguenza_3_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_NeoHookeanSiguenza.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of NeoHookeanSiguenza for dimension 3
         */
        template <typename T>
        class NeoHookeanSiguenza<T, 3> {
        public:
            static constexpr int Dim = 3;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "NeoHookeanSiguenza_3"; }

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    in.get("G", G);
                    in.get("K", K);
                }

                T G{1.0};
                T K{1.0};
            };

            UTOPIA_FUNCTION void hessian(const Params &params,
                                         const T *UTOPIA_RESTRICT f,
                                         const T *grad_test,
                                         const T *grad_trial,
                                         const T dx,
                                         T *UTOPIA_RESTRICT bf) const {
                const T G = params.G;
                const T K = params.K;

                using namespace utopia::device;
                // Automatically generated
                T x0 = f[4] * f[8];
                T x1 = f[5] * f[7];
                T x2 = x0 - x1;
                T x3 = f[5] * f[6];
                T x4 = f[3] * f[7];
                T x5 = f[3] * f[8];
                T x6 = f[4] * f[6];
                T x7 = f[0] * x0 - f[0] * x1 + f[1] * x3 - f[1] * x5 + f[2] * x4 - f[2] * x6;
                T x8 = K / pow(x7, 2);
                T x9 = x2 * x8;
                T x10 = log(x7);
                T x11 = -x10 * x2;
                T x12 = 2 * pow(x7, -0.66666666666666663);
                T x13 = 0.66666666666666663 * f[5] * f[7] - 0.66666666666666663 * x0;
                T x14 = pow(x7, -1.6666666666666665);
                T x15 = f[0] * x14;
                T x16 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) +
                        pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2);
                T x17 = x16 * pow(x7, -2.6666666666666665);
                T x18 = x17 * (1.6666666666666665 * f[5] * f[7] - 1.6666666666666665 * x0);
                T x19 = (1.0 / 2.0) * G;
                T x20 = -f[5] * f[6] + x5;
                T x21 = -x20;
                T x22 = x21 * x9;
                T x23 = x21 * x8;
                T x24 = -0.66666666666666663 * x3 + 0.66666666666666663 * x5;
                T x25 = 2 * x15;
                T x26 = f[1] * x14;
                T x27 = 2 * x13;
                T x28 = x24 * x25 + x26 * x27;
                T x29 = x4 - x6;
                T x30 = x29 * x9;
                T x31 = x29 * x8;
                T x32 = 0.66666666666666663 * f[4] * f[6] - 0.66666666666666663 * x4;
                T x33 = f[2] * x14;
                T x34 = x25 * x32 + x27 * x33;
                T x35 = x10 * x20;
                T x36 = x17 * (-1.6666666666666665 * x3 + 1.6666666666666665 * x5);
                T x37 = x23 * x29;
                T x38 = 2 * x26;
                T x39 = 2 * x24;
                T x40 = x32 * x38 + x33 * x39;
                T x41 = -x10 * x29;
                T x42 = x17 * (1.6666666666666665 * f[4] * f[6] - 1.6666666666666665 * x4);
                T x43 = f[1] * f[8];
                T x44 = -f[2] * f[7] + x43;
                T x45 = -x44;
                T x46 = x45 * x9;
                T x47 = x45 * x8;
                T x48 = f[2] * f[7];
                T x49 = 0.66666666666666663 * x43 - 0.66666666666666663 * x48;
                T x50 = x14 * x27;
                T x51 = f[3] * x50 + x25 * x49;
                T x52 = f[0] * f[8];
                T x53 = 0.66666666666666663 * f[2] * f[6] - 0.66666666666666663 * x52;
                T x54 = f[8] * x14;
                T x55 = 0.66666666666666663 * x16;
                T x56 = x54 * x55;
                T x57 = f[4] * x50 + x25 * x53 - x56;
                T x58 = -f[2] * f[6] + x52;
                T x59 = x58 * x8;
                T x60 = K * x10 / x7;
                T x61 = f[8] * x60;
                T x62 = x58 * x9 + x61;
                T x63 = f[0] * f[7];
                T x64 = f[1] * f[6];
                T x65 = 0.66666666666666663 * x63 - 0.66666666666666663 * x64;
                T x66 = x14 * x55;
                T x67 = f[7] * x66;
                T x68 = f[5] * x50 + x25 * x65 + x67;
                T x69 = -f[1] * f[6] + x63;
                T x70 = -x69;
                T x71 = x70 * x8;
                T x72 = f[7] * x60;
                T x73 = x70 * x9 - x72;
                T x74 = x23 * x58;
                T x75 = x14 * x39;
                T x76 = f[4] * x75 + x38 * x53;
                T x77 = f[6] * x66;
                T x78 = f[5] * x75 + x38 * x65 - x77;
                T x79 = f[6] * x60;
                T x80 = x23 * x70 + x79;
                T x81 = f[3] * x75 + x38 * x49 + x56;
                T x82 = x23 * x45 - x61;
                T x83 = x31 * x70;
                T x84 = 2 * x33;
                T x85 = 2 * x32;
                T x86 = x14 * x85;
                T x87 = f[5] * x86 + x65 * x84;
                T x88 = f[3] * x14;
                T x89 = x49 * x84 - x67 + x85 * x88;
                T x90 = x31 * x45 + x72;
                T x91 = f[4] * x86 + x53 * x84 + x77;
                T x92 = x31 * x58 - x79;
                T x93 = f[1] * f[5];
                T x94 = -f[2] * f[4] + x93;
                T x95 = x9 * x94;
                T x96 = x8 * x94;
                T x97 = 0.66666666666666663 * f[2] * f[4] - 0.66666666666666663 * x93;
                T x98 = f[6] * x50 + x25 * x97;
                T x99 = f[0] * f[4];
                T x100 = 0.66666666666666663 * f[1] * f[3] - 0.66666666666666663 * x99;
                T x101 = f[4] * x66;
                T x102 = f[8] * x50 + x100 * x25 - x101;
                T x103 = -f[1] * f[3] + x99;
                T x104 = x103 * x8;
                T x105 = f[4] * x60;
                T x106 = x103 * x9 + x105;
                T x107 = f[0] * f[5];
                T x108 = f[2] * f[3];
                T x109 = 0.66666666666666663 * x107 - 0.66666666666666663 * x108;
                T x110 = f[5] * x66;
                T x111 = f[7] * x50 + x109 * x25 + x110;
                T x112 = -f[2] * f[3] + x107;
                T x113 = -x112;
                T x114 = x113 * x8;
                T x115 = f[5] * x60;
                T x116 = x113 * x9 - x115;
                T x117 = x113 * x23;
                T x118 = f[7] * x75 + x109 * x38;
                T x119 = f[6] * x75 - x110 + x38 * x97;
                T x120 = x115 + x23 * x94;
                T x121 = f[3] * x66;
                T x122 = x100 * x38 + x121 + x39 * x54;
                T x123 = f[3] * x60;
                T x124 = x103 * x23 - x123;
                T x125 = x103 * x31;
                T x126 = x100 * x84 + x54 * x85;
                T x127 = f[7] * x86 + x109 * x84 - x121;
                T x128 = x113 * x31 + x123;
                T x129 = f[6] * x86 + x101 + x84 * x97;
                T x130 = -x105 + x31 * x94;
                T x131 = x10 * x44;
                T x132 = x17 * (1.6666666666666665 * x43 - 1.6666666666666665 * x48);
                T x133 = -x10 * x58;
                T x134 = x17 * (1.6666666666666665 * f[2] * f[6] - 1.6666666666666665 * x52);
                T x135 = x10 * x69;
                T x136 = x17 * (1.6666666666666665 * x63 - 1.6666666666666665 * x64);
                T x137 = x47 * x58;
                T x138 = 2 * x88;
                T x139 = f[4] * x14;
                T x140 = 2 * x49;
                T x141 = x138 * x53 + x139 * x140;
                T x142 = x47 * x70;
                T x143 = f[5] * x14;
                T x144 = x138 * x65 + x140 * x143;
                T x145 = x59 * x70;
                T x146 = 2 * x139;
                T x147 = 2 * x53;
                T x148 = x143 * x147 + x146 * x65;
                T x149 = x47 * x94;
                T x150 = x14 * x140;
                T x151 = f[6] * x150 + x138 * x97;
                T x152 = x33 * x55;
                T x153 = f[7] * x150 + x109 * x138 - x152;
                T x154 = f[2] * x60;
                T x155 = x113 * x47 + x154;
                T x156 = x26 * x55;
                T x157 = x100 * x138 + x140 * x54 + x156;
                T x158 = f[1] * x60;
                T x159 = x103 * x47 - x158;
                T x160 = x113 * x59;
                T x161 = x14 * x147;
                T x162 = f[7] * x161 + x109 * x146;
                T x163 = x15 * x55;
                T x164 = x100 * x146 + x147 * x54 - x163;
                T x165 = f[0] * x60;
                T x166 = x103 * x59 + x165;
                T x167 = f[6] * x161 + x146 * x97 + x152;
                T x168 = -x154 + x59 * x94;
                T x169 = x103 * x71;
                T x170 = 2 * x143;
                T x171 = 2 * x65;
                T x172 = x100 * x170 + x171 * x54;
                T x173 = f[6] * x14;
                T x174 = -x156 + x170 * x97 + x171 * x173;
                T x175 = x158 + x71 * x94;
                T x176 = f[7] * x14;
                T x177 = x109 * x170 + x163 + x171 * x176;
                T x178 = x113 * x71 - x165;
                T x179 = -x10 * x94;
                T x180 = x17 * (1.6666666666666665 * f[2] * f[4] - 1.6666666666666665 * x93);
                T x181 = x10 * x112;
                T x182 = x17 * (1.6666666666666665 * x107 - 1.6666666666666665 * x108);
                T x183 = -x10 * x103;
                T x184 = x17 * (1.6666666666666665 * f[1] * f[3] - 1.6666666666666665 * x99);
                T x185 = x113 * x96;
                T x186 = 2 * x173;
                T x187 = 2 * x97;
                T x188 = x109 * x186 + x176 * x187;
                T x189 = x103 * x96;
                T x190 = x100 * x186 + x187 * x54;
                T x191 = x103 * x114;
                T x192 = 2 * x100 * x176 + 2 * x109 * x54;
                bf[0] +=
                    dx *
                    (grad_test[0] *
                         (grad_trial[0] * (x11 * x9 + x19 * (x12 + 4 * x13 * x15 + x13 * x18) + pow(x2, 2) * x8) +
                          grad_trial[1] * (x11 * x23 + x19 * (x18 * x24 + x28) + x22) +
                          grad_trial[2] * (x11 * x31 + x19 * (x18 * x32 + x34) + x30)) +
                     grad_test[1] *
                         (grad_trial[0] * (x19 * (x13 * x36 + x28) + x22 + x35 * x9) +
                          grad_trial[1] * (x19 * (x12 + 4 * x24 * x26 + x24 * x36) + pow(x21, 2) * x8 + x23 * x35) +
                          grad_trial[2] * (x19 * (x32 * x36 + x40) + x31 * x35 + x37)) +
                     grad_test[2] *
                         (grad_trial[0] * (x19 * (x13 * x42 + x34) + x30 + x41 * x9) +
                          grad_trial[1] * (x19 * (x24 * x42 + x40) + x23 * x41 + x37) +
                          grad_trial[2] * (x19 * (x12 + 4 * x32 * x33 + x32 * x42) + pow(x29, 2) * x8 + x31 * x41)));
                bf[1] += dx * (grad_test[0] * (grad_trial[0] * (x11 * x47 + x19 * (x18 * x49 + x51) + x46) +
                                               grad_trial[1] * (x11 * x59 + x19 * (x18 * x53 + x57) + x62) +
                                               grad_trial[2] * (x11 * x71 + x19 * (x18 * x65 + x68) + x73)) +
                               grad_test[1] * (grad_trial[0] * (x19 * (x36 * x49 + x81) + x35 * x47 + x82) +
                                               grad_trial[1] * (x19 * (x36 * x53 + x76) + x35 * x59 + x74) +
                                               grad_trial[2] * (x19 * (x36 * x65 + x78) + x35 * x71 + x80)) +
                               grad_test[2] * (grad_trial[0] * (x19 * (x42 * x49 + x89) + x41 * x47 + x90) +
                                               grad_trial[1] * (x19 * (x42 * x53 + x91) + x41 * x59 + x92) +
                                               grad_trial[2] * (x19 * (x42 * x65 + x87) + x41 * x71 + x83)));
                bf[2] += dx * (grad_test[0] * (grad_trial[0] * (x11 * x96 + x19 * (x18 * x97 + x98) + x95) +
                                               grad_trial[1] * (x11 * x114 + x116 + x19 * (x109 * x18 + x111)) +
                                               grad_trial[2] * (x104 * x11 + x106 + x19 * (x100 * x18 + x102))) +
                               grad_test[1] * (grad_trial[0] * (x120 + x19 * (x119 + x36 * x97) + x35 * x96) +
                                               grad_trial[1] * (x114 * x35 + x117 + x19 * (x109 * x36 + x118)) +
                                               grad_trial[2] * (x104 * x35 + x124 + x19 * (x100 * x36 + x122))) +
                               grad_test[2] * (grad_trial[0] * (x130 + x19 * (x129 + x42 * x97) + x41 * x96) +
                                               grad_trial[1] * (x114 * x41 + x128 + x19 * (x109 * x42 + x127)) +
                                               grad_trial[2] * (x104 * x41 + x125 + x19 * (x100 * x42 + x126))));
                bf[3] += dx * (grad_test[0] * (grad_trial[0] * (x131 * x9 + x19 * (x13 * x132 + x51) + x46) +
                                               grad_trial[1] * (x131 * x23 + x19 * (x132 * x24 + x81) + x82) +
                                               grad_trial[2] * (x131 * x31 + x19 * (x132 * x32 + x89) + x90)) +
                               grad_test[1] * (grad_trial[0] * (x133 * x9 + x19 * (x13 * x134 + x57) + x62) +
                                               grad_trial[1] * (x133 * x23 + x19 * (x134 * x24 + x76) + x74) +
                                               grad_trial[2] * (x133 * x31 + x19 * (x134 * x32 + x91) + x92)) +
                               grad_test[2] * (grad_trial[0] * (x135 * x9 + x19 * (x13 * x136 + x68) + x73) +
                                               grad_trial[1] * (x135 * x23 + x19 * (x136 * x24 + x78) + x80) +
                                               grad_trial[2] * (x135 * x31 + x19 * (x136 * x32 + x87) + x83)));
                bf[4] +=
                    dx *
                    (grad_test[0] *
                         (grad_trial[0] * (x131 * x47 + x19 * (x12 + x132 * x49 + 4 * x49 * x88) + pow(x45, 2) * x8) +
                          grad_trial[1] * (x131 * x59 + x137 + x19 * (x132 * x53 + x141)) +
                          grad_trial[2] * (x131 * x71 + x142 + x19 * (x132 * x65 + x144))) +
                     grad_test[1] *
                         (grad_trial[0] * (x133 * x47 + x137 + x19 * (x134 * x49 + x141)) +
                          grad_trial[1] * (x133 * x59 + x19 * (x12 + x134 * x53 + 4 * x139 * x53) + pow(x58, 2) * x8) +
                          grad_trial[2] * (x133 * x71 + x145 + x19 * (x134 * x65 + x148))) +
                     grad_test[2] *
                         (grad_trial[0] * (x135 * x47 + x142 + x19 * (x136 * x49 + x144)) +
                          grad_trial[1] * (x135 * x59 + x145 + x19 * (x136 * x53 + x148)) +
                          grad_trial[2] * (x135 * x71 + x19 * (x12 + x136 * x65 + 4 * x143 * x65) + pow(x70, 2) * x8)));
                bf[5] += dx * (grad_test[0] * (grad_trial[0] * (x131 * x96 + x149 + x19 * (x132 * x97 + x151)) +
                                               grad_trial[1] * (x114 * x131 + x155 + x19 * (x109 * x132 + x153)) +
                                               grad_trial[2] * (x104 * x131 + x159 + x19 * (x100 * x132 + x157))) +
                               grad_test[1] * (grad_trial[0] * (x133 * x96 + x168 + x19 * (x134 * x97 + x167)) +
                                               grad_trial[1] * (x114 * x133 + x160 + x19 * (x109 * x134 + x162)) +
                                               grad_trial[2] * (x104 * x133 + x166 + x19 * (x100 * x134 + x164))) +
                               grad_test[2] * (grad_trial[0] * (x135 * x96 + x175 + x19 * (x136 * x97 + x174)) +
                                               grad_trial[1] * (x114 * x135 + x178 + x19 * (x109 * x136 + x177)) +
                                               grad_trial[2] * (x104 * x135 + x169 + x19 * (x100 * x136 + x172))));
                bf[6] += dx * (grad_test[0] * (grad_trial[0] * (x179 * x9 + x19 * (x13 * x180 + x98) + x95) +
                                               grad_trial[1] * (x120 + x179 * x23 + x19 * (x119 + x180 * x24)) +
                                               grad_trial[2] * (x130 + x179 * x31 + x19 * (x129 + x180 * x32))) +
                               grad_test[1] * (grad_trial[0] * (x116 + x181 * x9 + x19 * (x111 + x13 * x182)) +
                                               grad_trial[1] * (x117 + x181 * x23 + x19 * (x118 + x182 * x24)) +
                                               grad_trial[2] * (x128 + x181 * x31 + x19 * (x127 + x182 * x32))) +
                               grad_test[2] * (grad_trial[0] * (x106 + x183 * x9 + x19 * (x102 + x13 * x184)) +
                                               grad_trial[1] * (x124 + x183 * x23 + x19 * (x122 + x184 * x24)) +
                                               grad_trial[2] * (x125 + x183 * x31 + x19 * (x126 + x184 * x32))));
                bf[7] += dx * (grad_test[0] * (grad_trial[0] * (x149 + x179 * x47 + x19 * (x151 + x180 * x49)) +
                                               grad_trial[1] * (x168 + x179 * x59 + x19 * (x167 + x180 * x53)) +
                                               grad_trial[2] * (x175 + x179 * x71 + x19 * (x174 + x180 * x65))) +
                               grad_test[1] * (grad_trial[0] * (x155 + x181 * x47 + x19 * (x153 + x182 * x49)) +
                                               grad_trial[1] * (x160 + x181 * x59 + x19 * (x162 + x182 * x53)) +
                                               grad_trial[2] * (x178 + x181 * x71 + x19 * (x177 + x182 * x65))) +
                               grad_test[2] * (grad_trial[0] * (x159 + x183 * x47 + x19 * (x157 + x184 * x49)) +
                                               grad_trial[1] * (x166 + x183 * x59 + x19 * (x164 + x184 * x53)) +
                                               grad_trial[2] * (x169 + x183 * x71 + x19 * (x172 + x184 * x65))));
                bf[8] += dx * (grad_test[0] * (grad_trial[0] * (x179 * x96 + x19 * (x12 + 4 * x173 * x97 + x180 * x97) +
                                                                x8 * pow(x94, 2)) +
                                               grad_trial[1] * (x114 * x179 + x185 + x19 * (x109 * x180 + x188)) +
                                               grad_trial[2] * (x104 * x179 + x189 + x19 * (x100 * x180 + x190))) +
                               grad_test[1] * (grad_trial[0] * (x181 * x96 + x185 + x19 * (x182 * x97 + x188)) +
                                               grad_trial[1] * (pow(x113, 2) * x8 + x114 * x181 +
                                                                x19 * (4 * x109 * x176 + x109 * x182 + x12)) +
                                               grad_trial[2] * (x104 * x181 + x19 * (x100 * x182 + x192) + x191)) +
                               grad_test[2] * (grad_trial[0] * (x183 * x96 + x189 + x19 * (x184 * x97 + x190)) +
                                               grad_trial[1] * (x114 * x183 + x19 * (x109 * x184 + x192) + x191) +
                                               grad_trial[2] * (pow(x103, 2) * x8 + x104 * x183 +
                                                                x19 * (x100 * x184 + 4 * x100 * x54 + x12))));
            }

            UTOPIA_FUNCTION void gradient(const Params &params,
                                          const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                const T G = params.G;
                const T K = params.K;

                using namespace utopia::device;
                // Automatically generated
                T x0 = f[4] * f[8];
                T x1 = f[5] * f[7];
                T x2 = f[5] * f[6];
                T x3 = f[3] * f[7];
                T x4 = f[3] * f[8];
                T x5 = f[4] * f[6];
                T x6 = f[0] * x0 - f[0] * x1 + f[1] * x2 - f[1] * x4 + f[2] * x3 - f[2] * x5;
                T x7 = K * log(x6) / x6;
                T x8 = 2 * pow(x6, -0.66666666666666663);
                T x9 = pow(x6, -1.6666666666666665) *
                       (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) +
                        pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2));
                T x10 = (1.0 / 2.0) * G;
                T x11 = f[1] * f[8];
                T x12 = f[0] * f[8];
                T x13 = f[0] * f[7];
                T x14 = f[1] * f[5];
                T x15 = f[0] * f[5];
                T x16 = f[0] * f[4];
                lf[0] +=
                    dx *
                    (grad_test[0] *
                         (x10 * (f[0] * x8 + x9 * (0.66666666666666663 * f[5] * f[7] - 0.66666666666666663 * x0)) +
                          x7 * (x0 - x1)) +
                     grad_test[1] * (x10 * (f[1] * x8 + x9 * (-0.66666666666666663 * x2 + 0.66666666666666663 * x4)) +
                                     x7 * (f[5] * f[6] - x4)) +
                     grad_test[2] *
                         (x10 * (f[2] * x8 + x9 * (0.66666666666666663 * f[4] * f[6] - 0.66666666666666663 * x3)) +
                          x7 * (x3 - x5)));
                lf[1] +=
                    dx *
                    (grad_test[0] *
                         (x10 * (f[3] * x8 + x9 * (-0.66666666666666663 * f[2] * f[7] + 0.66666666666666663 * x11)) +
                          x7 * (f[2] * f[7] - x11)) +
                     grad_test[1] *
                         (x10 * (f[4] * x8 + x9 * (0.66666666666666663 * f[2] * f[6] - 0.66666666666666663 * x12)) +
                          x7 * (-f[2] * f[6] + x12)) +
                     grad_test[2] *
                         (x10 * (f[5] * x8 + x9 * (-0.66666666666666663 * f[1] * f[6] + 0.66666666666666663 * x13)) +
                          x7 * (f[1] * f[6] - x13)));
                lf[2] +=
                    dx *
                    (grad_test[0] *
                         (x10 * (f[6] * x8 + x9 * (0.66666666666666663 * f[2] * f[4] - 0.66666666666666663 * x14)) +
                          x7 * (-f[2] * f[4] + x14)) +
                     grad_test[1] *
                         (x10 * (f[7] * x8 + x9 * (-0.66666666666666663 * f[2] * f[3] + 0.66666666666666663 * x15)) +
                          x7 * (f[2] * f[3] - x15)) +
                     grad_test[2] *
                         (x10 * (f[8] * x8 + x9 * (0.66666666666666663 * f[1] * f[3] - 0.66666666666666663 * x16)) +
                          x7 * (-f[1] * f[3] + x16)));
            }

            UTOPIA_FUNCTION void value(const Params &params, const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                const T G = params.G;
                const T K = params.K;
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[0] * f[4] * f[8] - f[0] * f[5] * f[7] - f[1] * f[3] * f[8] + f[1] * f[5] * f[6] +
                       f[2] * f[3] * f[7] - f[2] * f[4] * f[6];
                e += dx * ((1.0 / 2.0) * G *
                               (pow(x0, -0.66666666666666663) *
                                    (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) +
                                     pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2)) -
                                3) +
                           (1.0 / 2.0) * K * pow(log(x0), 2));
            }

            UTOPIA_FUNCTION void eval(const Params &params,
                                      const T *UTOPIA_RESTRICT f,
                                      const T *grad_test,
                                      const T *grad_trial,
                                      const T dx,
                                      T &e,
                                      T *UTOPIA_RESTRICT lf,
                                      T *UTOPIA_RESTRICT bf) const {
                const T G = params.G;
                const T K = params.K;

                using namespace utopia::device;
                // Automatically generated
                T x0 = f[4] * f[8];
                T x1 = f[5] * f[6];
                T x2 = f[3] * f[7];
                T x3 = f[5] * f[7];
                T x4 = f[3] * f[8];
                T x5 = f[4] * f[6];
                T x6 = f[0] * x0 - f[0] * x3 + f[1] * x1 - f[1] * x4 + f[2] * x2 - f[2] * x5;
                T x7 = log(x6);
                T x8 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) +
                       pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2);
                T x9 = pow(x6, -0.66666666666666663);
                T x10 = (1.0 / 2.0) * G;
                T x11 = x0 - x3;
                T x12 = x11 * x7;
                T x13 = K / x6;
                T x14 = 2 * x9;
                T x15 = 0.66666666666666663 * f[5] * f[7] - 0.66666666666666663 * x0;
                T x16 = pow(x6, -1.6666666666666665);
                T x17 = x16 * x8;
                T x18 = -f[5] * f[6] + x4;
                T x19 = -x18;
                T x20 = x13 * x7;
                T x21 = -0.66666666666666663 * x1 + 0.66666666666666663 * x4;
                T x22 = x2 - x5;
                T x23 = 0.66666666666666663 * f[4] * f[6] - 0.66666666666666663 * x2;
                T x24 = f[1] * f[8];
                T x25 = -f[2] * f[7] + x24;
                T x26 = -x25;
                T x27 = f[2] * f[7];
                T x28 = 0.66666666666666663 * x24 - 0.66666666666666663 * x27;
                T x29 = f[0] * f[8];
                T x30 = -f[2] * f[6] + x29;
                T x31 = 0.66666666666666663 * f[2] * f[6] - 0.66666666666666663 * x29;
                T x32 = f[0] * f[7];
                T x33 = -f[1] * f[6] + x32;
                T x34 = -x33;
                T x35 = f[1] * f[6];
                T x36 = 0.66666666666666663 * x32 - 0.66666666666666663 * x35;
                T x37 = f[1] * f[5];
                T x38 = -f[2] * f[4] + x37;
                T x39 = 0.66666666666666663 * f[2] * f[4] - 0.66666666666666663 * x37;
                T x40 = f[0] * f[5];
                T x41 = -f[2] * f[3] + x40;
                T x42 = -x41;
                T x43 = f[2] * f[3];
                T x44 = 0.66666666666666663 * x40 - 0.66666666666666663 * x43;
                T x45 = f[0] * f[4];
                T x46 = -f[1] * f[3] + x45;
                T x47 = 0.66666666666666663 * f[1] * f[3] - 0.66666666666666663 * x45;
                T x48 = K / pow(x6, 2);
                T x49 = -x11;
                T x50 = x12 * x48;
                T x51 = f[0] * x16;
                T x52 = pow(x6, -2.6666666666666665) * x8;
                T x53 = x52 * (1.6666666666666665 * f[5] * f[7] - 1.6666666666666665 * x0);
                T x54 = x11 * x48;
                T x55 = x19 * x54;
                T x56 = x19 * x48;
                T x57 = x49 * x7;
                T x58 = 2 * x51;
                T x59 = f[1] * x16;
                T x60 = 2 * x15;
                T x61 = x21 * x58 + x59 * x60;
                T x62 = x22 * x54;
                T x63 = x22 * x48;
                T x64 = f[2] * x16;
                T x65 = x23 * x58 + x60 * x64;
                T x66 = x18 * x7;
                T x67 = x52 * (-1.6666666666666665 * x1 + 1.6666666666666665 * x4);
                T x68 = x22 * x56;
                T x69 = 2 * x59;
                T x70 = 2 * x21;
                T x71 = x23 * x69 + x64 * x70;
                T x72 = -x22;
                T x73 = x7 * x72;
                T x74 = x52 * (1.6666666666666665 * f[4] * f[6] - 1.6666666666666665 * x2);
                T x75 = x26 * x54;
                T x76 = x26 * x48;
                T x77 = x16 * x60;
                T x78 = f[3] * x77 + x28 * x58;
                T x79 = 0.66666666666666663 * x17;
                T x80 = f[8] * x79;
                T x81 = f[4] * x77 + x31 * x58 - x80;
                T x82 = x30 * x48;
                T x83 = f[8] * x20;
                T x84 = x30 * x54 + x83;
                T x85 = f[7] * x79;
                T x86 = f[5] * x77 + x36 * x58 + x85;
                T x87 = x34 * x48;
                T x88 = f[7] * x20;
                T x89 = x34 * x54 - x88;
                T x90 = x30 * x56;
                T x91 = x16 * x70;
                T x92 = f[4] * x91 + x31 * x69;
                T x93 = f[6] * x79;
                T x94 = f[5] * x91 + x36 * x69 - x93;
                T x95 = f[6] * x20;
                T x96 = x34 * x56 + x95;
                T x97 = f[3] * x91 + x28 * x69 + x80;
                T x98 = x26 * x56 - x83;
                T x99 = x34 * x63;
                T x100 = 2 * x64;
                T x101 = 2 * x23;
                T x102 = x101 * x16;
                T x103 = f[5] * x102 + x100 * x36;
                T x104 = f[3] * x16;
                T x105 = x100 * x28 + x101 * x104 - x85;
                T x106 = x26 * x63 + x88;
                T x107 = f[4] * x102 + x100 * x31 + x93;
                T x108 = x30 * x63 - x95;
                T x109 = x38 * x54;
                T x110 = x38 * x48;
                T x111 = f[6] * x77 + x39 * x58;
                T x112 = f[4] * x79;
                T x113 = f[8] * x77 - x112 + x47 * x58;
                T x114 = x46 * x48;
                T x115 = f[4] * x20;
                T x116 = x115 + x46 * x54;
                T x117 = f[5] * x79;
                T x118 = f[7] * x77 + x117 + x44 * x58;
                T x119 = x42 * x48;
                T x120 = f[5] * x20;
                T x121 = -x120 + x42 * x54;
                T x122 = x42 * x56;
                T x123 = f[7] * x91 + x44 * x69;
                T x124 = f[6] * x91 - x117 + x39 * x69;
                T x125 = x120 + x38 * x56;
                T x126 = f[3] * x79;
                T x127 = f[8] * x91 + x126 + x47 * x69;
                T x128 = f[3] * x20;
                T x129 = -x128 + x46 * x56;
                T x130 = x46 * x63;
                T x131 = f[8] * x102 + x100 * x47;
                T x132 = f[7] * x102 + x100 * x44 - x126;
                T x133 = x128 + x42 * x63;
                T x134 = f[6] * x102 + x100 * x39 + x112;
                T x135 = -x115 + x38 * x63;
                T x136 = x52 * (1.6666666666666665 * x24 - 1.6666666666666665 * x27);
                T x137 = x25 * x7;
                T x138 = -x30;
                T x139 = x138 * x7;
                T x140 = x52 * (1.6666666666666665 * f[2] * f[6] - 1.6666666666666665 * x29);
                T x141 = x33 * x7;
                T x142 = x52 * (1.6666666666666665 * x32 - 1.6666666666666665 * x35);
                T x143 = x30 * x76;
                T x144 = 2 * x104;
                T x145 = f[4] * x16;
                T x146 = 2 * x28;
                T x147 = x144 * x31 + x145 * x146;
                T x148 = x34 * x76;
                T x149 = f[5] * x16;
                T x150 = x144 * x36 + x146 * x149;
                T x151 = x34 * x82;
                T x152 = 2 * x145;
                T x153 = 2 * x31;
                T x154 = x149 * x153 + x152 * x36;
                T x155 = x38 * x76;
                T x156 = x146 * x16;
                T x157 = f[6] * x156 + x144 * x39;
                T x158 = f[2] * x79;
                T x159 = f[7] * x156 + x144 * x44 - x158;
                T x160 = f[2] * x20;
                T x161 = x160 + x42 * x76;
                T x162 = f[1] * x79;
                T x163 = f[8] * x156 + x144 * x47 + x162;
                T x164 = f[1] * x20;
                T x165 = -x164 + x46 * x76;
                T x166 = x42 * x82;
                T x167 = x153 * x16;
                T x168 = f[7] * x167 + x152 * x44;
                T x169 = f[0] * x79;
                T x170 = f[8] * x167 + x152 * x47 - x169;
                T x171 = f[0] * x20;
                T x172 = x171 + x46 * x82;
                T x173 = f[6] * x167 + x152 * x39 + x158;
                T x174 = -x160 + x38 * x82;
                T x175 = x46 * x87;
                T x176 = 2 * x149;
                T x177 = 2 * x36;
                T x178 = x16 * x177;
                T x179 = f[8] * x178 + x176 * x47;
                T x180 = f[6] * x16;
                T x181 = -x162 + x176 * x39 + x177 * x180;
                T x182 = x164 + x38 * x87;
                T x183 = f[7] * x178 + x169 + x176 * x44;
                T x184 = -x171 + x42 * x87;
                T x185 = -x38;
                T x186 = x52 * (1.6666666666666665 * f[2] * f[4] - 1.6666666666666665 * x37);
                T x187 = x185 * x7;
                T x188 = x41 * x7;
                T x189 = x52 * (1.6666666666666665 * x40 - 1.6666666666666665 * x43);
                T x190 = -x46;
                T x191 = x190 * x7;
                T x192 = x52 * (1.6666666666666665 * f[1] * f[3] - 1.6666666666666665 * x45);
                T x193 = x110 * x42;
                T x194 = 2 * x180;
                T x195 = f[7] * x16;
                T x196 = 2 * x39;
                T x197 = x194 * x44 + x195 * x196;
                T x198 = x110 * x46;
                T x199 = f[8] * x16;
                T x200 = x194 * x47 + x196 * x199;
                T x201 = x119 * x46;
                T x202 = 2 * x195 * x47 + 2 * x199 * x44;
                e += dx * ((1.0 / 2.0) * K * pow(x7, 2) + x10 * (x8 * x9 - 3));
                lf[0] += dx * (grad_test[0] * (x10 * (f[0] * x14 + x15 * x17) + x12 * x13) +
                               grad_test[1] * (x10 * (f[1] * x14 + x17 * x21) + x19 * x20) +
                               grad_test[2] * (x10 * (f[2] * x14 + x17 * x23) + x20 * x22));
                lf[1] += dx * (grad_test[0] * (x10 * (f[3] * x14 + x17 * x28) + x20 * x26) +
                               grad_test[1] * (x10 * (f[4] * x14 + x17 * x31) + x20 * x30) +
                               grad_test[2] * (x10 * (f[5] * x14 + x17 * x36) + x20 * x34));
                lf[2] += dx * (grad_test[0] * (x10 * (f[6] * x14 + x17 * x39) + x20 * x38) +
                               grad_test[1] * (x10 * (f[7] * x14 + x17 * x44) + x20 * x42) +
                               grad_test[2] * (x10 * (f[8] * x14 + x17 * x47) + x20 * x46));
                bf[0] +=
                    dx *
                    (grad_test[0] *
                         (grad_trial[0] * (x10 * (x14 + 4 * x15 * x51 + x15 * x53) + pow(x11, 2) * x48 + x49 * x50) +
                          grad_trial[1] * (x10 * (x21 * x53 + x61) + x55 + x56 * x57) +
                          grad_trial[2] * (x10 * (x23 * x53 + x65) + x57 * x63 + x62)) +
                     grad_test[1] *
                         (grad_trial[0] * (x10 * (x15 * x67 + x61) + x18 * x50 + x55) +
                          grad_trial[1] * (x10 * (x14 + 4 * x21 * x59 + x21 * x67) + pow(x19, 2) * x48 + x56 * x66) +
                          grad_trial[2] * (x10 * (x23 * x67 + x71) + x63 * x66 + x68)) +
                     grad_test[2] *
                         (grad_trial[0] * (x10 * (x15 * x74 + x65) + x50 * x72 + x62) +
                          grad_trial[1] * (x10 * (x21 * x74 + x71) + x56 * x73 + x68) +
                          grad_trial[2] * (x10 * (x14 + 4 * x23 * x64 + x23 * x74) + pow(x22, 2) * x48 + x63 * x73)));
                bf[1] += dx * (grad_test[0] * (grad_trial[0] * (x10 * (x28 * x53 + x78) + x57 * x76 + x75) +
                                               grad_trial[1] * (x10 * (x31 * x53 + x81) + x57 * x82 + x84) +
                                               grad_trial[2] * (x10 * (x36 * x53 + x86) + x57 * x87 + x89)) +
                               grad_test[1] * (grad_trial[0] * (x10 * (x28 * x67 + x97) + x66 * x76 + x98) +
                                               grad_trial[1] * (x10 * (x31 * x67 + x92) + x66 * x82 + x90) +
                                               grad_trial[2] * (x10 * (x36 * x67 + x94) + x66 * x87 + x96)) +
                               grad_test[2] * (grad_trial[0] * (x10 * (x105 + x28 * x74) + x106 + x73 * x76) +
                                               grad_trial[1] * (x10 * (x107 + x31 * x74) + x108 + x73 * x82) +
                                               grad_trial[2] * (x10 * (x103 + x36 * x74) + x73 * x87 + x99)));
                bf[2] += dx * (grad_test[0] * (grad_trial[0] * (x10 * (x111 + x39 * x53) + x109 + x110 * x57) +
                                               grad_trial[1] * (x10 * (x118 + x44 * x53) + x119 * x57 + x121) +
                                               grad_trial[2] * (x10 * (x113 + x47 * x53) + x114 * x57 + x116)) +
                               grad_test[1] * (grad_trial[0] * (x10 * (x124 + x39 * x67) + x110 * x66 + x125) +
                                               grad_trial[1] * (x10 * (x123 + x44 * x67) + x119 * x66 + x122) +
                                               grad_trial[2] * (x10 * (x127 + x47 * x67) + x114 * x66 + x129)) +
                               grad_test[2] * (grad_trial[0] * (x10 * (x134 + x39 * x74) + x110 * x73 + x135) +
                                               grad_trial[1] * (x10 * (x132 + x44 * x74) + x119 * x73 + x133) +
                                               grad_trial[2] * (x10 * (x131 + x47 * x74) + x114 * x73 + x130)));
                bf[3] += dx * (grad_test[0] * (grad_trial[0] * (x10 * (x136 * x15 + x78) + x25 * x50 + x75) +
                                               grad_trial[1] * (x10 * (x136 * x21 + x97) + x137 * x56 + x98) +
                                               grad_trial[2] * (x10 * (x105 + x136 * x23) + x106 + x137 * x63)) +
                               grad_test[1] * (grad_trial[0] * (x10 * (x140 * x15 + x81) + x138 * x50 + x84) +
                                               grad_trial[1] * (x10 * (x140 * x21 + x92) + x139 * x56 + x90) +
                                               grad_trial[2] * (x10 * (x107 + x140 * x23) + x108 + x139 * x63)) +
                               grad_test[2] * (grad_trial[0] * (x10 * (x142 * x15 + x86) + x33 * x50 + x89) +
                                               grad_trial[1] * (x10 * (x142 * x21 + x94) + x141 * x56 + x96) +
                                               grad_trial[2] * (x10 * (x103 + x142 * x23) + x141 * x63 + x99)));
                bf[4] +=
                    dx *
                    (grad_test[0] *
                         (grad_trial[0] * (x10 * (4 * x104 * x28 + x136 * x28 + x14) + x137 * x76 + pow(x26, 2) * x48) +
                          grad_trial[1] * (x10 * (x136 * x31 + x147) + x137 * x82 + x143) +
                          grad_trial[2] * (x10 * (x136 * x36 + x150) + x137 * x87 + x148)) +
                     grad_test[1] *
                         (grad_trial[0] * (x10 * (x140 * x28 + x147) + x139 * x76 + x143) +
                          grad_trial[1] * (x10 * (x14 + x140 * x31 + 4 * x145 * x31) + x139 * x82 + pow(x30, 2) * x48) +
                          grad_trial[2] * (x10 * (x140 * x36 + x154) + x139 * x87 + x151)) +
                     grad_test[2] * (grad_trial[0] * (x10 * (x142 * x28 + x150) + x141 * x76 + x148) +
                                     grad_trial[1] * (x10 * (x142 * x31 + x154) + x141 * x82 + x151) +
                                     grad_trial[2] *
                                         (x10 * (x14 + x142 * x36 + 4 * x149 * x36) + x141 * x87 + pow(x34, 2) * x48)));
                bf[5] += dx * (grad_test[0] * (grad_trial[0] * (x10 * (x136 * x39 + x157) + x110 * x137 + x155) +
                                               grad_trial[1] * (x10 * (x136 * x44 + x159) + x119 * x137 + x161) +
                                               grad_trial[2] * (x10 * (x136 * x47 + x163) + x114 * x137 + x165)) +
                               grad_test[1] * (grad_trial[0] * (x10 * (x140 * x39 + x173) + x110 * x139 + x174) +
                                               grad_trial[1] * (x10 * (x140 * x44 + x168) + x119 * x139 + x166) +
                                               grad_trial[2] * (x10 * (x140 * x47 + x170) + x114 * x139 + x172)) +
                               grad_test[2] * (grad_trial[0] * (x10 * (x142 * x39 + x181) + x110 * x141 + x182) +
                                               grad_trial[1] * (x10 * (x142 * x44 + x183) + x119 * x141 + x184) +
                                               grad_trial[2] * (x10 * (x142 * x47 + x179) + x114 * x141 + x175)));
                bf[6] += dx * (grad_test[0] * (grad_trial[0] * (x10 * (x111 + x15 * x186) + x109 + x185 * x50) +
                                               grad_trial[1] * (x10 * (x124 + x186 * x21) + x125 + x187 * x56) +
                                               grad_trial[2] * (x10 * (x134 + x186 * x23) + x135 + x187 * x63)) +
                               grad_test[1] * (grad_trial[0] * (x10 * (x118 + x15 * x189) + x121 + x41 * x50) +
                                               grad_trial[1] * (x10 * (x123 + x189 * x21) + x122 + x188 * x56) +
                                               grad_trial[2] * (x10 * (x132 + x189 * x23) + x133 + x188 * x63)) +
                               grad_test[2] * (grad_trial[0] * (x10 * (x113 + x15 * x192) + x116 + x190 * x50) +
                                               grad_trial[1] * (x10 * (x127 + x192 * x21) + x129 + x191 * x56) +
                                               grad_trial[2] * (x10 * (x131 + x192 * x23) + x130 + x191 * x63)));
                bf[7] += dx * (grad_test[0] * (grad_trial[0] * (x10 * (x157 + x186 * x28) + x155 + x187 * x76) +
                                               grad_trial[1] * (x10 * (x173 + x186 * x31) + x174 + x187 * x82) +
                                               grad_trial[2] * (x10 * (x181 + x186 * x36) + x182 + x187 * x87)) +
                               grad_test[1] * (grad_trial[0] * (x10 * (x159 + x189 * x28) + x161 + x188 * x76) +
                                               grad_trial[1] * (x10 * (x168 + x189 * x31) + x166 + x188 * x82) +
                                               grad_trial[2] * (x10 * (x183 + x189 * x36) + x184 + x188 * x87)) +
                               grad_test[2] * (grad_trial[0] * (x10 * (x163 + x192 * x28) + x165 + x191 * x76) +
                                               grad_trial[1] * (x10 * (x170 + x192 * x31) + x172 + x191 * x82) +
                                               grad_trial[2] * (x10 * (x179 + x192 * x36) + x175 + x191 * x87)));
                bf[8] += dx * (grad_test[0] * (grad_trial[0] * (x10 * (x14 + 4 * x180 * x39 + x186 * x39) +
                                                                x110 * x187 + pow(x38, 2) * x48) +
                                               grad_trial[1] * (x10 * (x186 * x44 + x197) + x119 * x187 + x193) +
                                               grad_trial[2] * (x10 * (x186 * x47 + x200) + x114 * x187 + x198)) +
                               grad_test[1] * (grad_trial[0] * (x10 * (x189 * x39 + x197) + x110 * x188 + x193) +
                                               grad_trial[1] * (x10 * (x14 + x189 * x44 + 4 * x195 * x44) +
                                                                x119 * x188 + pow(x42, 2) * x48) +
                                               grad_trial[2] * (x10 * (x189 * x47 + x202) + x114 * x188 + x201)) +
                               grad_test[2] * (grad_trial[0] * (x10 * (x192 * x39 + x200) + x110 * x191 + x198) +
                                               grad_trial[1] * (x10 * (x192 * x44 + x202) + x119 * x191 + x201) +
                                               grad_trial[2] * (x10 * (x14 + x192 * x47 + 4 * x199 * x47) +
                                                                x114 * x191 + pow(x46, 2) * x48)));
            }

            UTOPIA_FUNCTION void apply(const Params &params,
                                       const T *UTOPIA_RESTRICT f,
                                       const T *grad_test,
                                       const T *disp_grad,
                                       const T dx,
                                       T *UTOPIA_RESTRICT res) const {
                const T G = params.G;
                const T K = params.K;

                using namespace utopia::device;
                // Automatically generated
                T x0 = f[4] * f[8];
                T x1 = f[5] * f[7];
                T x2 = x0 - x1;
                T x3 = f[5] * f[6];
                T x4 = f[3] * f[7];
                T x5 = f[3] * f[8];
                T x6 = f[4] * f[6];
                T x7 = f[0] * x0 - f[0] * x1 + f[1] * x3 - f[1] * x5 + f[2] * x4 - f[2] * x6;
                T x8 = K / pow(x7, 2);
                T x9 = x2 * x8;
                T x10 = log(x7);
                T x11 = -x10 * x2;
                T x12 = 2 * pow(x7, -0.66666666666666663);
                T x13 = 0.66666666666666663 * f[5] * f[7] - 0.66666666666666663 * x0;
                T x14 = pow(x7, -1.6666666666666665);
                T x15 = f[0] * x14;
                T x16 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) +
                        pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2);
                T x17 = x16 * pow(x7, -2.6666666666666665);
                T x18 = x17 * (1.6666666666666665 * f[5] * f[7] - 1.6666666666666665 * x0);
                T x19 = (1.0 / 2.0) * G;
                T x20 = -f[5] * f[6] + x5;
                T x21 = -x20;
                T x22 = x21 * x9;
                T x23 = x21 * x8;
                T x24 = -0.66666666666666663 * x3 + 0.66666666666666663 * x5;
                T x25 = 2 * x15;
                T x26 = f[1] * x14;
                T x27 = 2 * x13;
                T x28 = x24 * x25 + x26 * x27;
                T x29 = x4 - x6;
                T x30 = x29 * x9;
                T x31 = x29 * x8;
                T x32 = 0.66666666666666663 * f[4] * f[6] - 0.66666666666666663 * x4;
                T x33 = x14 * x27;
                T x34 = f[2] * x33 + x25 * x32;
                T x35 = f[1] * f[8];
                T x36 = -f[2] * f[7] + x35;
                T x37 = -x36;
                T x38 = x37 * x9;
                T x39 = x37 * x8;
                T x40 = f[2] * f[7];
                T x41 = 0.66666666666666663 * x35 - 0.66666666666666663 * x40;
                T x42 = f[3] * x33 + x25 * x41;
                T x43 = f[1] * f[5];
                T x44 = -f[2] * f[4] + x43;
                T x45 = x44 * x9;
                T x46 = x44 * x8;
                T x47 = 0.66666666666666663 * f[2] * f[4] - 0.66666666666666663 * x43;
                T x48 = f[6] * x33 + x25 * x47;
                T x49 = f[0] * f[8];
                T x50 = 0.66666666666666663 * f[2] * f[6] - 0.66666666666666663 * x49;
                T x51 = f[8] * x14;
                T x52 = 0.66666666666666663 * x16;
                T x53 = x51 * x52;
                T x54 = f[4] * x33 + x25 * x50 - x53;
                T x55 = -f[2] * f[6] + x49;
                T x56 = x55 * x8;
                T x57 = K * x10 / x7;
                T x58 = f[8] * x57;
                T x59 = x55 * x9 + x58;
                T x60 = f[0] * f[4];
                T x61 = 0.66666666666666663 * f[1] * f[3] - 0.66666666666666663 * x60;
                T x62 = x14 * x52;
                T x63 = f[4] * x62;
                T x64 = f[8] * x33 + x25 * x61 - x63;
                T x65 = -f[1] * f[3] + x60;
                T x66 = x65 * x8;
                T x67 = f[4] * x57;
                T x68 = x65 * x9 + x67;
                T x69 = f[0] * f[7];
                T x70 = f[1] * f[6];
                T x71 = 0.66666666666666663 * x69 - 0.66666666666666663 * x70;
                T x72 = f[7] * x62;
                T x73 = f[5] * x33 + x25 * x71 + x72;
                T x74 = -f[1] * f[6] + x69;
                T x75 = -x74;
                T x76 = x75 * x8;
                T x77 = f[7] * x57;
                T x78 = x75 * x9 - x77;
                T x79 = f[0] * f[5];
                T x80 = f[2] * f[3];
                T x81 = 0.66666666666666663 * x79 - 0.66666666666666663 * x80;
                T x82 = f[5] * x62;
                T x83 = f[7] * x33 + x25 * x81 + x82;
                T x84 = -f[2] * f[3] + x79;
                T x85 = -x84;
                T x86 = x8 * x85;
                T x87 = f[5] * x57;
                T x88 = x85 * x9 - x87;
                T x89 = x10 * x20;
                T x90 = x17 * (-1.6666666666666665 * x3 + 1.6666666666666665 * x5);
                T x91 = x23 * x29;
                T x92 = 2 * x26;
                T x93 = 2 * x24;
                T x94 = f[2] * x14;
                T x95 = x32 * x92 + x93 * x94;
                T x96 = x23 * x55;
                T x97 = x14 * x93;
                T x98 = f[4] * x97 + x50 * x92;
                T x99 = x23 * x85;
                T x100 = f[7] * x97 + x81 * x92;
                T x101 = f[6] * x62;
                T x102 = f[5] * x97 - x101 + x71 * x92;
                T x103 = f[6] * x57;
                T x104 = x103 + x23 * x75;
                T x105 = f[6] * x97 + x47 * x92 - x82;
                T x106 = x23 * x44 + x87;
                T x107 = f[3] * x97 + x41 * x92 + x53;
                T x108 = x23 * x37 - x58;
                T x109 = f[3] * x62;
                T x110 = x109 + x51 * x93 + x61 * x92;
                T x111 = f[3] * x57;
                T x112 = -x111 + x23 * x65;
                T x113 = -x10 * x29;
                T x114 = x17 * (1.6666666666666665 * f[4] * f[6] - 1.6666666666666665 * x4);
                T x115 = x31 * x75;
                T x116 = 2 * x94;
                T x117 = 2 * x32;
                T x118 = x117 * x14;
                T x119 = f[5] * x118 + x116 * x71;
                T x120 = x31 * x65;
                T x121 = x116 * x61 + x117 * x51;
                T x122 = f[3] * x14;
                T x123 = x116 * x41 + x117 * x122 - x72;
                T x124 = x31 * x37 + x77;
                T x125 = f[7] * x118 - x109 + x116 * x81;
                T x126 = x111 + x31 * x85;
                T x127 = f[4] * x118 + x101 + x116 * x50;
                T x128 = -x103 + x31 * x55;
                T x129 = f[6] * x118 + x116 * x47 + x63;
                T x130 = x31 * x44 - x67;
                T x131 = x10 * x36;
                T x132 = x17 * (1.6666666666666665 * x35 - 1.6666666666666665 * x40);
                T x133 = x39 * x55;
                T x134 = 2 * x122;
                T x135 = 2 * x41;
                T x136 = f[4] * x14;
                T x137 = x134 * x50 + x135 * x136;
                T x138 = x39 * x75;
                T x139 = x135 * x14;
                T x140 = f[5] * x139 + x134 * x71;
                T x141 = x39 * x44;
                T x142 = f[6] * x139 + x134 * x47;
                T x143 = f[2] * x62;
                T x144 = f[7] * x139 + x134 * x81 - x143;
                T x145 = f[2] * x57;
                T x146 = x145 + x39 * x85;
                T x147 = x26 * x52;
                T x148 = x134 * x61 + x135 * x51 + x147;
                T x149 = f[1] * x57;
                T x150 = -x149 + x39 * x65;
                T x151 = -x10 * x55;
                T x152 = x17 * (1.6666666666666665 * f[2] * f[6] - 1.6666666666666665 * x49);
                T x153 = x56 * x75;
                T x154 = 2 * x136;
                T x155 = 2 * x50;
                T x156 = f[5] * x14;
                T x157 = x154 * x71 + x155 * x156;
                T x158 = x56 * x85;
                T x159 = x14 * x155;
                T x160 = f[7] * x159 + x154 * x81;
                T x161 = x15 * x52;
                T x162 = x154 * x61 + x155 * x51 - x161;
                T x163 = f[0] * x57;
                T x164 = x163 + x56 * x65;
                T x165 = f[6] * x159 + x143 + x154 * x47;
                T x166 = -x145 + x44 * x56;
                T x167 = x10 * x74;
                T x168 = x17 * (1.6666666666666665 * x69 - 1.6666666666666665 * x70);
                T x169 = x65 * x76;
                T x170 = 2 * x156;
                T x171 = 2 * x71;
                T x172 = x170 * x61 + x171 * x51;
                T x173 = f[6] * x14;
                T x174 = -x147 + x170 * x47 + x171 * x173;
                T x175 = x149 + x44 * x76;
                T x176 = f[7] * x14;
                T x177 = x161 + x170 * x81 + x171 * x176;
                T x178 = -x163 + x76 * x85;
                T x179 = -x10 * x44;
                T x180 = x17 * (1.6666666666666665 * f[2] * f[4] - 1.6666666666666665 * x43);
                T x181 = x46 * x85;
                T x182 = 2 * x173;
                T x183 = 2 * x47;
                T x184 = x176 * x183 + x182 * x81;
                T x185 = x46 * x65;
                T x186 = x182 * x61 + x183 * x51;
                T x187 = x10 * x84;
                T x188 = x17 * (1.6666666666666665 * x79 - 1.6666666666666665 * x80);
                T x189 = x65 * x86;
                T x190 = 2 * x176 * x61 + 2 * x51 * x81;
                T x191 = -x10 * x65;
                T x192 = x17 * (1.6666666666666665 * f[1] * f[3] - 1.6666666666666665 * x60);
                res[0] +=
                    dx * (grad_test[0] *
                              (disp_grad[0] * (x11 * x9 + x19 * (x12 + 4 * x13 * x15 + x13 * x18) + pow(x2, 2) * x8) +
                               disp_grad[1] * (x11 * x23 + x19 * (x18 * x24 + x28) + x22) +
                               disp_grad[2] * (x11 * x31 + x19 * (x18 * x32 + x34) + x30) +
                               disp_grad[3] * (x11 * x39 + x19 * (x18 * x41 + x42) + x38) +
                               disp_grad[4] * (x11 * x56 + x19 * (x18 * x50 + x54) + x59) +
                               disp_grad[5] * (x11 * x76 + x19 * (x18 * x71 + x73) + x78) +
                               disp_grad[6] * (x11 * x46 + x19 * (x18 * x47 + x48) + x45) +
                               disp_grad[7] * (x11 * x86 + x19 * (x18 * x81 + x83) + x88) +
                               disp_grad[8] * (x11 * x66 + x19 * (x18 * x61 + x64) + x68)) +
                          grad_test[1] *
                              (disp_grad[0] * (x19 * (x13 * x90 + x28) + x22 + x89 * x9) +
                               disp_grad[1] * (x19 * (x12 + 4 * x24 * x26 + x24 * x90) + pow(x21, 2) * x8 + x23 * x89) +
                               disp_grad[2] * (x19 * (x32 * x90 + x95) + x31 * x89 + x91) +
                               disp_grad[3] * (x108 + x19 * (x107 + x41 * x90) + x39 * x89) +
                               disp_grad[4] * (x19 * (x50 * x90 + x98) + x56 * x89 + x96) +
                               disp_grad[5] * (x104 + x19 * (x102 + x71 * x90) + x76 * x89) +
                               disp_grad[6] * (x106 + x19 * (x105 + x47 * x90) + x46 * x89) +
                               disp_grad[7] * (x19 * (x100 + x81 * x90) + x86 * x89 + x99) +
                               disp_grad[8] * (x112 + x19 * (x110 + x61 * x90) + x66 * x89)) +
                          grad_test[2] * (disp_grad[0] * (x113 * x9 + x19 * (x114 * x13 + x34) + x30) +
                                          disp_grad[1] * (x113 * x23 + x19 * (x114 * x24 + x95) + x91) +
                                          disp_grad[2] * (x113 * x31 + x19 * (x114 * x32 + x12 + 4 * x32 * x94) +
                                                          pow(x29, 2) * x8) +
                                          disp_grad[3] * (x113 * x39 + x124 + x19 * (x114 * x41 + x123)) +
                                          disp_grad[4] * (x113 * x56 + x128 + x19 * (x114 * x50 + x127)) +
                                          disp_grad[5] * (x113 * x76 + x115 + x19 * (x114 * x71 + x119)) +
                                          disp_grad[6] * (x113 * x46 + x130 + x19 * (x114 * x47 + x129)) +
                                          disp_grad[7] * (x113 * x86 + x126 + x19 * (x114 * x81 + x125)) +
                                          disp_grad[8] * (x113 * x66 + x120 + x19 * (x114 * x61 + x121))));
                res[1] +=
                    dx *
                    (grad_test[0] *
                         (disp_grad[0] * (x131 * x9 + x19 * (x13 * x132 + x42) + x38) +
                          disp_grad[1] * (x108 + x131 * x23 + x19 * (x107 + x132 * x24)) +
                          disp_grad[2] * (x124 + x131 * x31 + x19 * (x123 + x132 * x32)) +
                          disp_grad[3] * (x131 * x39 + x19 * (x12 + 4 * x122 * x41 + x132 * x41) + pow(x37, 2) * x8) +
                          disp_grad[4] * (x131 * x56 + x133 + x19 * (x132 * x50 + x137)) +
                          disp_grad[5] * (x131 * x76 + x138 + x19 * (x132 * x71 + x140)) +
                          disp_grad[6] * (x131 * x46 + x141 + x19 * (x132 * x47 + x142)) +
                          disp_grad[7] * (x131 * x86 + x146 + x19 * (x132 * x81 + x144)) +
                          disp_grad[8] * (x131 * x66 + x150 + x19 * (x132 * x61 + x148))) +
                     grad_test[1] *
                         (disp_grad[0] * (x151 * x9 + x19 * (x13 * x152 + x54) + x59) +
                          disp_grad[1] * (x151 * x23 + x19 * (x152 * x24 + x98) + x96) +
                          disp_grad[2] * (x128 + x151 * x31 + x19 * (x127 + x152 * x32)) +
                          disp_grad[3] * (x133 + x151 * x39 + x19 * (x137 + x152 * x41)) +
                          disp_grad[4] * (x151 * x56 + x19 * (x12 + 4 * x136 * x50 + x152 * x50) + pow(x55, 2) * x8) +
                          disp_grad[5] * (x151 * x76 + x153 + x19 * (x152 * x71 + x157)) +
                          disp_grad[6] * (x151 * x46 + x166 + x19 * (x152 * x47 + x165)) +
                          disp_grad[7] * (x151 * x86 + x158 + x19 * (x152 * x81 + x160)) +
                          disp_grad[8] * (x151 * x66 + x164 + x19 * (x152 * x61 + x162))) +
                     grad_test[2] *
                         (disp_grad[0] * (x167 * x9 + x19 * (x13 * x168 + x73) + x78) +
                          disp_grad[1] * (x104 + x167 * x23 + x19 * (x102 + x168 * x24)) +
                          disp_grad[2] * (x115 + x167 * x31 + x19 * (x119 + x168 * x32)) +
                          disp_grad[3] * (x138 + x167 * x39 + x19 * (x140 + x168 * x41)) +
                          disp_grad[4] * (x153 + x167 * x56 + x19 * (x157 + x168 * x50)) +
                          disp_grad[5] * (x167 * x76 + x19 * (x12 + 4 * x156 * x71 + x168 * x71) + pow(x75, 2) * x8) +
                          disp_grad[6] * (x167 * x46 + x175 + x19 * (x168 * x47 + x174)) +
                          disp_grad[7] * (x167 * x86 + x178 + x19 * (x168 * x81 + x177)) +
                          disp_grad[8] * (x167 * x66 + x169 + x19 * (x168 * x61 + x172))));
                res[2] +=
                    dx *
                    (grad_test[0] *
                         (disp_grad[0] * (x179 * x9 + x19 * (x13 * x180 + x48) + x45) +
                          disp_grad[1] * (x106 + x179 * x23 + x19 * (x105 + x180 * x24)) +
                          disp_grad[2] * (x130 + x179 * x31 + x19 * (x129 + x180 * x32)) +
                          disp_grad[3] * (x141 + x179 * x39 + x19 * (x142 + x180 * x41)) +
                          disp_grad[4] * (x166 + x179 * x56 + x19 * (x165 + x180 * x50)) +
                          disp_grad[5] * (x175 + x179 * x76 + x19 * (x174 + x180 * x71)) +
                          disp_grad[6] * (x179 * x46 + x19 * (x12 + 4 * x173 * x47 + x180 * x47) + pow(x44, 2) * x8) +
                          disp_grad[7] * (x179 * x86 + x181 + x19 * (x180 * x81 + x184)) +
                          disp_grad[8] * (x179 * x66 + x185 + x19 * (x180 * x61 + x186))) +
                     grad_test[1] *
                         (disp_grad[0] * (x187 * x9 + x19 * (x13 * x188 + x83) + x88) +
                          disp_grad[1] * (x187 * x23 + x19 * (x100 + x188 * x24) + x99) +
                          disp_grad[2] * (x126 + x187 * x31 + x19 * (x125 + x188 * x32)) +
                          disp_grad[3] * (x146 + x187 * x39 + x19 * (x144 + x188 * x41)) +
                          disp_grad[4] * (x158 + x187 * x56 + x19 * (x160 + x188 * x50)) +
                          disp_grad[5] * (x178 + x187 * x76 + x19 * (x177 + x188 * x71)) +
                          disp_grad[6] * (x181 + x187 * x46 + x19 * (x184 + x188 * x47)) +
                          disp_grad[7] * (x187 * x86 + x19 * (x12 + 4 * x176 * x81 + x188 * x81) + x8 * pow(x85, 2)) +
                          disp_grad[8] * (x187 * x66 + x189 + x19 * (x188 * x61 + x190))) +
                     grad_test[2] *
                         (disp_grad[0] * (x19 * (x13 * x192 + x64) + x191 * x9 + x68) +
                          disp_grad[1] * (x112 + x19 * (x110 + x192 * x24) + x191 * x23) +
                          disp_grad[2] * (x120 + x19 * (x121 + x192 * x32) + x191 * x31) +
                          disp_grad[3] * (x150 + x19 * (x148 + x192 * x41) + x191 * x39) +
                          disp_grad[4] * (x164 + x19 * (x162 + x192 * x50) + x191 * x56) +
                          disp_grad[5] * (x169 + x19 * (x172 + x192 * x71) + x191 * x76) +
                          disp_grad[6] * (x185 + x19 * (x186 + x192 * x47) + x191 * x46) +
                          disp_grad[7] * (x189 + x19 * (x190 + x192 * x81) + x191 * x86) +
                          disp_grad[8] * (x19 * (x12 + x192 * x61 + 4 * x51 * x61) + x191 * x66 + pow(x65, 2) * x8)));
            }
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_NeoHookeanSiguenza_3_IMPL_hpp
