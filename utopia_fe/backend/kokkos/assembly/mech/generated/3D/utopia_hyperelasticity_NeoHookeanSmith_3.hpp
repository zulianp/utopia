#ifndef UTOPIA_TPL_HYPERELASTICITY_NeoHookeanSmith_3_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeoHookeanSmith_3_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_NeoHookeanSmith.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of NeoHookeanSmith for dimension 3
         */
        template <typename T>
        class NeoHookeanSmith<T, 3> {
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

            NeoHookeanSmith(const Params &params) {
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
                T x0 = f[3] * f[8];
                T x1 = f[5] * f[6];
                T x2 = -2 * x0 + 2 * x1;
                T x3 = f[4] * f[8];
                T x4 = f[5] * f[7];
                T x5 = (11.0 / 12.0) * lambda;
                T x6 = x5 * (x3 - x4);
                T x7 = pow(f[0], 2);
                T x8 = pow(f[1], 2);
                T x9 = pow(f[2], 2);
                T x10 = pow(f[3], 2);
                T x11 = pow(f[4], 2);
                T x12 = pow(f[5], 2);
                T x13 = pow(f[6], 2);
                T x14 = pow(f[7], 2);
                T x15 = pow(f[8], 2);
                T x16 = x10 + x11 + x12 + x13 + x14 + x15 + x7 + x8 + x9 + 1;
                T x17 = (8.0 / 3.0) * mu / pow(x16, 2);
                T x18 = f[0] * x17;
                T x19 = f[1] * x18;
                T x20 = f[3] * f[7];
                T x21 = f[4] * f[6];
                T x22 = 2 * x20 - 2 * x21;
                T x23 = f[2] * x18;
                T x24 = 2 * x3 - 2 * x4;
                T x25 = (4.0 / 3.0) * mu;
                T x26 = x25 - x25 / x16;
                T x27 = x5 * (-x0 + x1);
                T x28 = f[1] * x17;
                T x29 = f[2] * x28;
                T x30 = x5 * (x20 - x21);
                T x31 = f[2] * f[7];
                T x32 = f[1] * f[8];
                T x33 = x5 * (x31 - x32);
                T x34 = f[3] * x18;
                T x35 = (11.0 / 6.0) * lambda *
                        (f[0] * x3 - f[0] * x4 - f[1] * x0 + f[1] * x1 + f[2] * x20 - f[2] * x21 - 1 -
                         6.0 / 11.0 * mu / lambda);
                T x36 = f[8] * x35;
                T x37 = f[1] * f[3];
                T x38 = x17 * x37 - x36;
                T x39 = f[7] * x35;
                T x40 = f[2] * f[3];
                T x41 = x17 * x40 + x39;
                T x42 = f[0] * f[8];
                T x43 = f[2] * f[6];
                T x44 = x5 * (x42 - x43);
                T x45 = f[4] * x28;
                T x46 = f[0] * f[4];
                T x47 = x17 * x46 + x36;
                T x48 = f[6] * x35;
                T x49 = f[2] * f[4];
                T x50 = x17 * x49 - x48;
                T x51 = f[1] * f[6];
                T x52 = f[0] * f[7];
                T x53 = x5 * (x51 - x52);
                T x54 = f[2] * x17;
                T x55 = f[5] * x54;
                T x56 = f[0] * f[5];
                T x57 = x17 * x56 - x39;
                T x58 = f[1] * f[5];
                T x59 = x17 * x58 + x48;
                T x60 = x5 * (-x49 + x58);
                T x61 = f[6] * x18;
                T x62 = f[5] * x35;
                T x63 = x17 * x51 + x62;
                T x64 = f[4] * x35;
                T x65 = x17 * x43 - x64;
                T x66 = x5 * (x40 - x56);
                T x67 = f[7] * x28;
                T x68 = x17 * x52 - x62;
                T x69 = f[3] * x35;
                T x70 = x17 * x31 + x69;
                T x71 = x5 * (-x37 + x46);
                T x72 = f[8] * x54;
                T x73 = x17 * x42 + x64;
                T x74 = x17 * x32 - x69;
                T x75 = 2 * x31 - 2 * x32;
                T x76 = 2 * x42 - 2 * x43;
                T x77 = 2 * x51 - 2 * x52;
                T x78 = f[3] * x17;
                T x79 = f[4] * x78;
                T x80 = f[5] * x78;
                T x81 = f[4] * x17;
                T x82 = f[5] * x81;
                T x83 = f[6] * x78;
                T x84 = f[2] * x35;
                T x85 = x17 * x21 - x84;
                T x86 = f[1] * x35;
                T x87 = x1 * x17 + x86;
                T x88 = f[7] * x81;
                T x89 = x17 * x20 + x84;
                T x90 = f[0] * x35;
                T x91 = x17 * x4 - x90;
                T x92 = f[8] * x17;
                T x93 = f[5] * x92;
                T x94 = x0 * x17 - x86;
                T x95 = x17 * x3 + x90;
                T x96 = -2 * x49 + 2 * x58;
                T x97 = 2 * x40 - 2 * x56;
                T x98 = -2 * x37 + 2 * x46;
                T x99 = f[6] * f[7] * x17;
                T x100 = f[6] * x92;
                T x101 = f[7] * x92;
                bf[0] += dx * (grad_trial[0] * (grad_test[0] * (x17 * x7 + x24 * x6 + x26) +
                                                grad_test[1] * (x19 + x2 * x6) + grad_test[2] * (x22 * x6 + x23)) +
                               grad_trial[1] *
                                   (grad_test[0] * (x19 + x24 * x27) + grad_test[1] * (x17 * x8 + x2 * x27 + x26) +
                                    grad_test[2] * (x22 * x27 + x29)) +
                               grad_trial[2] * (grad_test[0] * (x23 + x24 * x30) + grad_test[1] * (x2 * x30 + x29) +
                                                grad_test[2] * (x17 * x9 + x22 * x30 + x26)));
                bf[1] += dx * (grad_trial[0] * (grad_test[0] * (x24 * x33 + x34) + grad_test[1] * (x2 * x33 + x38) +
                                                grad_test[2] * (x22 * x33 + x41)) +
                               grad_trial[1] * (grad_test[0] * (x24 * x44 + x47) + grad_test[1] * (x2 * x44 + x45) +
                                                grad_test[2] * (x22 * x44 + x50)) +
                               grad_trial[2] * (grad_test[0] * (x24 * x53 + x57) + grad_test[1] * (x2 * x53 + x59) +
                                                grad_test[2] * (x22 * x53 + x55)));
                bf[2] += dx * (grad_trial[0] * (grad_test[0] * (x24 * x60 + x61) + grad_test[1] * (x2 * x60 + x63) +
                                                grad_test[2] * (x22 * x60 + x65)) +
                               grad_trial[1] * (grad_test[0] * (x24 * x66 + x68) + grad_test[1] * (x2 * x66 + x67) +
                                                grad_test[2] * (x22 * x66 + x70)) +
                               grad_trial[2] * (grad_test[0] * (x24 * x71 + x73) + grad_test[1] * (x2 * x71 + x74) +
                                                grad_test[2] * (x22 * x71 + x72)));
                bf[3] += dx * (grad_trial[0] * (grad_test[0] * (x34 + x6 * x75) + grad_test[1] * (x47 + x6 * x76) +
                                                grad_test[2] * (x57 + x6 * x77)) +
                               grad_trial[1] * (grad_test[0] * (x27 * x75 + x38) + grad_test[1] * (x27 * x76 + x45) +
                                                grad_test[2] * (x27 * x77 + x59)) +
                               grad_trial[2] * (grad_test[0] * (x30 * x75 + x41) + grad_test[1] * (x30 * x76 + x50) +
                                                grad_test[2] * (x30 * x77 + x55)));
                bf[4] += dx * (grad_trial[0] * (grad_test[0] * (x10 * x17 + x26 + x33 * x75) +
                                                grad_test[1] * (x33 * x76 + x79) + grad_test[2] * (x33 * x77 + x80)) +
                               grad_trial[1] *
                                   (grad_test[0] * (x44 * x75 + x79) + grad_test[1] * (x11 * x17 + x26 + x44 * x76) +
                                    grad_test[2] * (x44 * x77 + x82)) +
                               grad_trial[2] * (grad_test[0] * (x53 * x75 + x80) + grad_test[1] * (x53 * x76 + x82) +
                                                grad_test[2] * (x12 * x17 + x26 + x53 * x77)));
                bf[5] += dx * (grad_trial[0] * (grad_test[0] * (x60 * x75 + x83) + grad_test[1] * (x60 * x76 + x85) +
                                                grad_test[2] * (x60 * x77 + x87)) +
                               grad_trial[1] * (grad_test[0] * (x66 * x75 + x89) + grad_test[1] * (x66 * x76 + x88) +
                                                grad_test[2] * (x66 * x77 + x91)) +
                               grad_trial[2] * (grad_test[0] * (x71 * x75 + x94) + grad_test[1] * (x71 * x76 + x95) +
                                                grad_test[2] * (x71 * x77 + x93)));
                bf[6] += dx * (grad_trial[0] * (grad_test[0] * (x6 * x96 + x61) + grad_test[1] * (x6 * x97 + x68) +
                                                grad_test[2] * (x6 * x98 + x73)) +
                               grad_trial[1] * (grad_test[0] * (x27 * x96 + x63) + grad_test[1] * (x27 * x97 + x67) +
                                                grad_test[2] * (x27 * x98 + x74)) +
                               grad_trial[2] * (grad_test[0] * (x30 * x96 + x65) + grad_test[1] * (x30 * x97 + x70) +
                                                grad_test[2] * (x30 * x98 + x72)));
                bf[7] += dx * (grad_trial[0] * (grad_test[0] * (x33 * x96 + x83) + grad_test[1] * (x33 * x97 + x89) +
                                                grad_test[2] * (x33 * x98 + x94)) +
                               grad_trial[1] * (grad_test[0] * (x44 * x96 + x85) + grad_test[1] * (x44 * x97 + x88) +
                                                grad_test[2] * (x44 * x98 + x95)) +
                               grad_trial[2] * (grad_test[0] * (x53 * x96 + x87) + grad_test[1] * (x53 * x97 + x91) +
                                                grad_test[2] * (x53 * x98 + x93)));
                bf[8] += dx * (grad_trial[0] * (grad_test[0] * (x13 * x17 + x26 + x60 * x96) +
                                                grad_test[1] * (x60 * x97 + x99) + grad_test[2] * (x100 + x60 * x98)) +
                               grad_trial[1] *
                                   (grad_test[0] * (x66 * x96 + x99) + grad_test[1] * (x14 * x17 + x26 + x66 * x97) +
                                    grad_test[2] * (x101 + x66 * x98)) +
                               grad_trial[2] * (grad_test[0] * (x100 + x71 * x96) + grad_test[1] * (x101 + x71 * x97) +
                                                grad_test[2] * (x15 * x17 + x26 + x71 * x98)));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = (4.0 / 3.0) * mu;
                T x1 = f[0] * x0;
                T x2 = 1.0 / (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) +
                              pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2) + 1);
                T x3 = f[4] * f[8];
                T x4 = f[5] * f[7];
                T x5 = f[5] * f[6];
                T x6 = f[3] * f[7];
                T x7 = f[3] * f[8];
                T x8 = f[4] * f[6];
                T x9 = (11.0 / 12.0) * lambda *
                       (f[0] * x3 - f[0] * x4 + f[1] * x5 - f[1] * x7 + f[2] * x6 - f[2] * x8 - 1 -
                        6.0 / 11.0 * mu / lambda);
                T x10 = f[1] * x0;
                T x11 = f[2] * x0;
                T x12 = f[3] * x0;
                T x13 = 2 * f[8];
                T x14 = 2 * f[2];
                T x15 = f[4] * x0;
                T x16 = f[5] * x0;
                T x17 = 2 * f[0];
                T x18 = 2 * f[1];
                T x19 = f[6] * x0;
                T x20 = f[7] * x0;
                T x21 = f[8] * x0;
                lf[0] += dx * (grad_test[0] * (-x1 * x2 + x1 + x9 * (2 * x3 - 2 * x4)) +
                               grad_test[1] * (-x10 * x2 + x10 + x9 * (2 * x5 - 2 * x7)) +
                               grad_test[2] * (-x11 * x2 + x11 + x9 * (2 * x6 - 2 * x8)));
                lf[1] += dx * (grad_test[0] * (-x12 * x2 + x12 + x9 * (-f[1] * x13 + f[7] * x14)) +
                               grad_test[1] * (-x15 * x2 + x15 + x9 * (f[0] * x13 - f[6] * x14)) +
                               grad_test[2] * (-x16 * x2 + x16 + x9 * (f[6] * x18 - f[7] * x17)));
                lf[2] += dx * (grad_test[0] * (-x19 * x2 + x19 + x9 * (-f[4] * x14 + f[5] * x18)) +
                               grad_test[1] * (-x2 * x20 + x20 + x9 * (f[3] * x14 - f[5] * x17)) +
                               grad_test[2] * (-x2 * x21 + x21 + x9 * (-f[3] * x18 + f[4] * x17)));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) +
                       pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2);
                T x1 = (2.0 / 3.0) * mu;
                e += dx * ((11.0 / 12.0) * lambda *
                               pow(f[0] * f[4] * f[8] - f[0] * f[5] * f[7] - f[1] * f[3] * f[8] + f[1] * f[5] * f[6] +
                                       f[2] * f[3] * f[7] - f[2] * f[4] * f[6] - 1 - 6.0 / 11.0 * mu / lambda,
                                   2) +
                           x1 * (x0 - 3) - x1 * log(x0 + 1));
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
                T x0 = pow(f[0], 2);
                T x1 = pow(f[1], 2);
                T x2 = pow(f[2], 2);
                T x3 = pow(f[3], 2);
                T x4 = pow(f[4], 2);
                T x5 = pow(f[5], 2);
                T x6 = pow(f[6], 2);
                T x7 = pow(f[7], 2);
                T x8 = pow(f[8], 2);
                T x9 = x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8;
                T x10 = (2.0 / 3.0) * mu;
                T x11 = x9 + 1;
                T x12 = f[4] * f[8];
                T x13 = f[5] * f[6];
                T x14 = f[3] * f[7];
                T x15 = f[5] * f[7];
                T x16 = f[3] * f[8];
                T x17 = f[4] * f[6];
                T x18 = f[0] * x12 - f[0] * x15 + f[1] * x13 - f[1] * x16 + f[2] * x14 - f[2] * x17 - 1 -
                        6.0 / 11.0 * mu / lambda;
                T x19 = (11.0 / 12.0) * lambda;
                T x20 = (4.0 / 3.0) * mu;
                T x21 = f[0] * x20;
                T x22 = 1.0 / x11;
                T x23 = 2 * x12 - 2 * x15;
                T x24 = x18 * x19;
                T x25 = f[1] * x20;
                T x26 = 2 * x13 - 2 * x16;
                T x27 = f[2] * x20;
                T x28 = 2 * x14 - 2 * x17;
                T x29 = f[3] * x20;
                T x30 = f[1] * f[8];
                T x31 = f[2] * f[7];
                T x32 = -2 * x30 + 2 * x31;
                T x33 = f[4] * x20;
                T x34 = f[0] * f[8];
                T x35 = f[2] * f[6];
                T x36 = 2 * x34 - 2 * x35;
                T x37 = f[5] * x20;
                T x38 = f[0] * f[7];
                T x39 = f[1] * f[6];
                T x40 = -2 * x38 + 2 * x39;
                T x41 = f[6] * x20;
                T x42 = f[1] * f[5];
                T x43 = f[2] * f[4];
                T x44 = 2 * x42 - 2 * x43;
                T x45 = f[7] * x20;
                T x46 = f[0] * f[5];
                T x47 = f[2] * f[3];
                T x48 = -2 * x46 + 2 * x47;
                T x49 = f[8] * x20;
                T x50 = f[0] * f[4];
                T x51 = f[1] * f[3];
                T x52 = 2 * x50 - 2 * x51;
                T x53 = x19 * (x12 - x15);
                T x54 = (8.0 / 3.0) * mu / pow(x11, 2);
                T x55 = f[0] * x54;
                T x56 = f[1] * x55;
                T x57 = f[2] * x55;
                T x58 = -x20 * x22 + x20;
                T x59 = x19 * (x13 - x16);
                T x60 = f[1] * x54;
                T x61 = f[2] * x60;
                T x62 = x19 * (x14 - x17);
                T x63 = x19 * (-x30 + x31);
                T x64 = f[3] * x55;
                T x65 = (11.0 / 6.0) * lambda * x18;
                T x66 = f[8] * x65;
                T x67 = x51 * x54 - x66;
                T x68 = f[7] * x65;
                T x69 = x47 * x54 + x68;
                T x70 = x19 * (x34 - x35);
                T x71 = f[4] * x60;
                T x72 = x50 * x54 + x66;
                T x73 = f[6] * x65;
                T x74 = x43 * x54 - x73;
                T x75 = x19 * (-x38 + x39);
                T x76 = f[2] * x54;
                T x77 = f[5] * x76;
                T x78 = x46 * x54 - x68;
                T x79 = x42 * x54 + x73;
                T x80 = x19 * (x42 - x43);
                T x81 = f[6] * x55;
                T x82 = f[5] * x65;
                T x83 = x39 * x54 + x82;
                T x84 = f[4] * x65;
                T x85 = x35 * x54 - x84;
                T x86 = x19 * (-x46 + x47);
                T x87 = f[7] * x60;
                T x88 = x38 * x54 - x82;
                T x89 = f[3] * x65;
                T x90 = x31 * x54 + x89;
                T x91 = x19 * (x50 - x51);
                T x92 = f[8] * x76;
                T x93 = x34 * x54 + x84;
                T x94 = x30 * x54 - x89;
                T x95 = f[3] * x54;
                T x96 = f[4] * x95;
                T x97 = f[5] * x95;
                T x98 = f[4] * x54;
                T x99 = f[5] * x98;
                T x100 = f[6] * x95;
                T x101 = f[2] * x65;
                T x102 = -x101 + x17 * x54;
                T x103 = f[1] * x65;
                T x104 = x103 + x13 * x54;
                T x105 = f[7] * x98;
                T x106 = x101 + x14 * x54;
                T x107 = f[0] * x65;
                T x108 = -x107 + x15 * x54;
                T x109 = f[8] * x54;
                T x110 = f[5] * x109;
                T x111 = -x103 + x16 * x54;
                T x112 = x107 + x12 * x54;
                T x113 = f[6] * f[7] * x54;
                T x114 = f[6] * x109;
                T x115 = f[7] * x109;
                e += dx * (x10 * (x9 - 3) - x10 * log(x11) + pow(x18, 2) * x19);
                lf[0] += dx * (grad_test[0] * (-x21 * x22 + x21 + x23 * x24) +
                               grad_test[1] * (-x22 * x25 + x24 * x26 + x25) +
                               grad_test[2] * (-x22 * x27 + x24 * x28 + x27));
                lf[1] += dx * (grad_test[0] * (-x22 * x29 + x24 * x32 + x29) +
                               grad_test[1] * (-x22 * x33 + x24 * x36 + x33) +
                               grad_test[2] * (-x22 * x37 + x24 * x40 + x37));
                lf[2] += dx * (grad_test[0] * (-x22 * x41 + x24 * x44 + x41) +
                               grad_test[1] * (-x22 * x45 + x24 * x48 + x45) +
                               grad_test[2] * (-x22 * x49 + x24 * x52 + x49));
                bf[0] += dx * (grad_trial[0] * (grad_test[0] * (x0 * x54 + x23 * x53 + x58) +
                                                grad_test[1] * (x26 * x53 + x56) + grad_test[2] * (x28 * x53 + x57)) +
                               grad_trial[1] *
                                   (grad_test[0] * (x23 * x59 + x56) + grad_test[1] * (x1 * x54 + x26 * x59 + x58) +
                                    grad_test[2] * (x28 * x59 + x61)) +
                               grad_trial[2] * (grad_test[0] * (x23 * x62 + x57) + grad_test[1] * (x26 * x62 + x61) +
                                                grad_test[2] * (x2 * x54 + x28 * x62 + x58)));
                bf[1] += dx * (grad_trial[0] * (grad_test[0] * (x23 * x63 + x64) + grad_test[1] * (x26 * x63 + x67) +
                                                grad_test[2] * (x28 * x63 + x69)) +
                               grad_trial[1] * (grad_test[0] * (x23 * x70 + x72) + grad_test[1] * (x26 * x70 + x71) +
                                                grad_test[2] * (x28 * x70 + x74)) +
                               grad_trial[2] * (grad_test[0] * (x23 * x75 + x78) + grad_test[1] * (x26 * x75 + x79) +
                                                grad_test[2] * (x28 * x75 + x77)));
                bf[2] += dx * (grad_trial[0] * (grad_test[0] * (x23 * x80 + x81) + grad_test[1] * (x26 * x80 + x83) +
                                                grad_test[2] * (x28 * x80 + x85)) +
                               grad_trial[1] * (grad_test[0] * (x23 * x86 + x88) + grad_test[1] * (x26 * x86 + x87) +
                                                grad_test[2] * (x28 * x86 + x90)) +
                               grad_trial[2] * (grad_test[0] * (x23 * x91 + x93) + grad_test[1] * (x26 * x91 + x94) +
                                                grad_test[2] * (x28 * x91 + x92)));
                bf[3] += dx * (grad_trial[0] * (grad_test[0] * (x32 * x53 + x64) + grad_test[1] * (x36 * x53 + x72) +
                                                grad_test[2] * (x40 * x53 + x78)) +
                               grad_trial[1] * (grad_test[0] * (x32 * x59 + x67) + grad_test[1] * (x36 * x59 + x71) +
                                                grad_test[2] * (x40 * x59 + x79)) +
                               grad_trial[2] * (grad_test[0] * (x32 * x62 + x69) + grad_test[1] * (x36 * x62 + x74) +
                                                grad_test[2] * (x40 * x62 + x77)));
                bf[4] += dx * (grad_trial[0] * (grad_test[0] * (x3 * x54 + x32 * x63 + x58) +
                                                grad_test[1] * (x36 * x63 + x96) + grad_test[2] * (x40 * x63 + x97)) +
                               grad_trial[1] *
                                   (grad_test[0] * (x32 * x70 + x96) + grad_test[1] * (x36 * x70 + x4 * x54 + x58) +
                                    grad_test[2] * (x40 * x70 + x99)) +
                               grad_trial[2] * (grad_test[0] * (x32 * x75 + x97) + grad_test[1] * (x36 * x75 + x99) +
                                                grad_test[2] * (x40 * x75 + x5 * x54 + x58)));
                bf[5] += dx * (grad_trial[0] * (grad_test[0] * (x100 + x32 * x80) + grad_test[1] * (x102 + x36 * x80) +
                                                grad_test[2] * (x104 + x40 * x80)) +
                               grad_trial[1] * (grad_test[0] * (x106 + x32 * x86) + grad_test[1] * (x105 + x36 * x86) +
                                                grad_test[2] * (x108 + x40 * x86)) +
                               grad_trial[2] * (grad_test[0] * (x111 + x32 * x91) + grad_test[1] * (x112 + x36 * x91) +
                                                grad_test[2] * (x110 + x40 * x91)));
                bf[6] += dx * (grad_trial[0] * (grad_test[0] * (x44 * x53 + x81) + grad_test[1] * (x48 * x53 + x88) +
                                                grad_test[2] * (x52 * x53 + x93)) +
                               grad_trial[1] * (grad_test[0] * (x44 * x59 + x83) + grad_test[1] * (x48 * x59 + x87) +
                                                grad_test[2] * (x52 * x59 + x94)) +
                               grad_trial[2] * (grad_test[0] * (x44 * x62 + x85) + grad_test[1] * (x48 * x62 + x90) +
                                                grad_test[2] * (x52 * x62 + x92)));
                bf[7] += dx * (grad_trial[0] * (grad_test[0] * (x100 + x44 * x63) + grad_test[1] * (x106 + x48 * x63) +
                                                grad_test[2] * (x111 + x52 * x63)) +
                               grad_trial[1] * (grad_test[0] * (x102 + x44 * x70) + grad_test[1] * (x105 + x48 * x70) +
                                                grad_test[2] * (x112 + x52 * x70)) +
                               grad_trial[2] * (grad_test[0] * (x104 + x44 * x75) + grad_test[1] * (x108 + x48 * x75) +
                                                grad_test[2] * (x110 + x52 * x75)));
                bf[8] += dx * (grad_trial[0] * (grad_test[0] * (x44 * x80 + x54 * x6 + x58) +
                                                grad_test[1] * (x113 + x48 * x80) + grad_test[2] * (x114 + x52 * x80)) +
                               grad_trial[1] *
                                   (grad_test[0] * (x113 + x44 * x86) + grad_test[1] * (x48 * x86 + x54 * x7 + x58) +
                                    grad_test[2] * (x115 + x52 * x86)) +
                               grad_trial[2] * (grad_test[0] * (x114 + x44 * x91) + grad_test[1] * (x115 + x48 * x91) +
                                                grad_test[2] * (x52 * x91 + x54 * x8 + x58)));
            }

            T mu{1.0};
            T lambda{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_NeoHookeanSmith_3_IMPL_hpp
