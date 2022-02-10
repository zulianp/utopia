#ifndef UTOPIA_TPL_HYPERELASTICITY_NeoHookeanOgden_3_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeoHookeanOgden_3_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_NeoHookeanOgden.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of NeoHookeanOgden for dimension 3
         */
        template <typename T>
        class NeoHookeanOgden<T, 3> {
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

            NeoHookeanOgden(const Params &params) {
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
                T x7 = log(x6);
                T x8 = lambda * x7;
                T x9 = mu - x8;
                T x10 = lambda + x9;
                T x11 = x0 - x3;
                T x12 = x10 * x11;
                T x13 = x2 - x5;
                T x14 = grad_test[2] * x13;
                T x15 = -x1 + x4;
                T x16 = x12 * x15;
                T x17 = pow(x11, 2);
                T x18 = lambda * x17;
                T x19 = pow(x6, 2);
                T x20 = mu * x19;
                T x21 = x10 * x15;
                T x22 = pow(x13, 2);
                T x23 = pow(x15, 2);
                T x24 = dx / x19;
                T x25 = f[1] * f[8] - f[2] * f[7];
                T x26 = x10 * x25;
                T x27 = grad_test[0] * x11;
                T x28 = -x26 * x27;
                T x29 = x6 * x9;
                T x30 = f[8] * x29 + x15 * x26;
                T x31 = f[7] * x29 + x13 * x26;
                T x32 = f[0] * f[8] - f[2] * f[6];
                T x33 = x10 * x32;
                T x34 = grad_test[1] * x15;
                T x35 = -x33 * x34;
                T x36 = f[6] * x29 + x13 * x33;
                T x37 = x6 * (-mu + x8);
                T x38 = f[8] * x37 + x11 * x33;
                T x39 = f[0] * f[7] - f[1] * f[6];
                T x40 = x10 * x39;
                T x41 = -x14 * x40;
                T x42 = f[6] * x37 + x15 * x40;
                T x43 = f[7] * x37 + x11 * x40;
                T x44 = f[1] * f[5] - f[2] * f[4];
                T x45 = x10 * x44;
                T x46 = x27 * x45;
                T x47 = f[4] * x29 + x13 * x45;
                T x48 = f[5] * x29 + x15 * x45;
                T x49 = f[0] * f[5] - f[2] * f[3];
                T x50 = x10 * x49;
                T x51 = x34 * x50;
                T x52 = f[3] * x29 + x13 * x50;
                T x53 = f[5] * x37 + x11 * x50;
                T x54 = f[0] * f[4] - f[1] * f[3];
                T x55 = x10 * x54;
                T x56 = x14 * x55;
                T x57 = f[4] * x37 + x11 * x55;
                T x58 = f[3] * x37 + x15 * x55;
                T x59 = grad_test[2] * x39;
                T x60 = x26 * x32;
                T x61 = pow(x25, 2);
                T x62 = grad_test[0] * x26;
                T x63 = grad_test[1] * x33;
                T x64 = pow(x39, 2);
                T x65 = pow(x32, 2);
                T x66 = -x44 * x62;
                T x67 = f[2] * x29 + x33 * x44;
                T x68 = f[1] * x29 + x40 * x44;
                T x69 = -x49 * x63;
                T x70 = f[0] * x29 + x40 * x49;
                T x71 = f[2] * x37 + x26 * x49;
                T x72 = grad_test[2] * x54;
                T x73 = -x40 * x72;
                T x74 = f[0] * x37 + x33 * x54;
                T x75 = f[1] * x37 + x26 * x54;
                T x76 = x45 * x49;
                T x77 = pow(x44, 2);
                T x78 = pow(x54, 2);
                T x79 = pow(x49, 2);
                bf[0] +=
                    x24 * (grad_trial[0] *
                               (grad_test[0] * (mu * x17 - x18 * x7 + x18 + x20) - grad_test[1] * x16 + x12 * x14) +
                           grad_trial[1] * (-grad_test[0] * x16 +
                                            grad_test[1] * (lambda * x23 + mu * x23 + x20 - x23 * x8) - x14 * x21) +
                           grad_trial[2] * (grad_test[0] * x12 * x13 - grad_test[1] * x13 * x21 +
                                            grad_test[2] * (lambda * x22 + mu * x22 + x20 - x22 * x8)));
                bf[1] += x24 * (grad_trial[0] * (grad_test[1] * x30 - grad_test[2] * x31 + x28) +
                                grad_trial[1] * (grad_test[0] * x38 + grad_test[2] * x36 + x35) +
                                grad_trial[2] * (-grad_test[0] * x43 + grad_test[1] * x42 + x41));
                bf[2] += x24 * (grad_trial[0] * (-grad_test[1] * x48 + grad_test[2] * x47 + x46) +
                                grad_trial[1] * (-grad_test[0] * x53 - grad_test[2] * x52 + x51) +
                                grad_trial[2] * (grad_test[0] * x57 - grad_test[1] * x58 + x56));
                bf[3] += x24 * (grad_trial[0] * (grad_test[1] * x38 - grad_test[2] * x43 + x28) +
                                grad_trial[1] * (grad_test[0] * x30 + grad_test[2] * x42 + x35) +
                                grad_trial[2] * (-grad_test[0] * x31 + grad_test[1] * x36 + x41));
                bf[4] +=
                    x24 * (grad_trial[0] * (grad_test[0] * (lambda * x61 + mu * x61 + x20 - x61 * x8) -
                                            grad_test[1] * x60 + x26 * x59) +
                           grad_trial[1] * (-grad_test[0] * x60 +
                                            grad_test[1] * (lambda * x65 + mu * x65 + x20 - x65 * x8) - x33 * x59) +
                           grad_trial[2] *
                               (grad_test[2] * (lambda * x64 + mu * x64 + x20 - x64 * x8) + x39 * x62 - x39 * x63));
                bf[5] += x24 * (grad_trial[0] * (grad_test[1] * x67 - grad_test[2] * x68 + x66) +
                                grad_trial[1] * (grad_test[0] * x71 + grad_test[2] * x70 + x69) +
                                grad_trial[2] * (-grad_test[0] * x75 + grad_test[1] * x74 + x73));
                bf[6] += x24 * (grad_trial[0] * (-grad_test[1] * x53 + grad_test[2] * x57 + x46) +
                                grad_trial[1] * (-grad_test[0] * x48 - grad_test[2] * x58 + x51) +
                                grad_trial[2] * (grad_test[0] * x47 - grad_test[1] * x52 + x56));
                bf[7] += x24 * (grad_trial[0] * (grad_test[1] * x71 - grad_test[2] * x75 + x66) +
                                grad_trial[1] * (grad_test[0] * x67 + grad_test[2] * x74 + x69) +
                                grad_trial[2] * (-grad_test[0] * x68 + grad_test[1] * x70 + x73));
                bf[8] +=
                    x24 * (grad_trial[0] * (grad_test[0] * (lambda * x77 + mu * x77 + x20 - x77 * x8) -
                                            grad_test[1] * x76 + x45 * x72) +
                           grad_trial[1] * (-grad_test[0] * x76 +
                                            grad_test[1] * (lambda * x79 + mu * x79 + x20 - x79 * x8) - x50 * x72) +
                           grad_trial[2] * (grad_test[0] * x45 * x54 - grad_test[1] * x50 * x54 +
                                            grad_test[2] * (lambda * x78 + mu * x78 + x20 - x78 * x8)));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
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
                T x8 = mu * x7;
                T x9 = lambda * log(x7);
                T x10 = -x3 + x5;
                T x11 = x4 - x6;
                T x12 = dx / x7;
                T x13 = f[1] * f[8] - f[2] * f[7];
                T x14 = f[0] * f[8] - f[2] * f[6];
                T x15 = f[0] * f[7] - f[1] * f[6];
                T x16 = f[1] * f[5] - f[2] * f[4];
                T x17 = f[0] * f[5] - f[2] * f[3];
                T x18 = f[0] * f[4] - f[1] * f[3];
                lf[0] += x12 * (grad_test[0] * (f[0] * x8 - mu * x2 + x2 * x9) +
                                grad_test[1] * (f[1] * x8 + mu * x10 - x10 * x9) +
                                grad_test[2] * (f[2] * x8 - mu * x11 + x11 * x9));
                lf[1] += x12 * (grad_test[0] * (f[3] * x8 + mu * x13 - x13 * x9) +
                                grad_test[1] * (f[4] * x8 - mu * x14 + x14 * x9) +
                                grad_test[2] * (f[5] * x8 + mu * x15 - x15 * x9));
                lf[2] += x12 * (grad_test[0] * (f[6] * x8 - mu * x16 + x16 * x9) +
                                grad_test[1] * (f[7] * x8 + mu * x17 - x17 * x9) +
                                grad_test[2] * (f[8] * x8 - mu * x18 + x18 * x9));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = log(f[0] * f[4] * f[8] - f[0] * f[5] * f[7] - f[1] * f[3] * f[8] + f[1] * f[5] * f[6] +
                           f[2] * f[3] * f[7] - f[2] * f[4] * f[6]);
                e += dx * ((1.0 / 2.0) * lambda * pow(x0, 2) - mu * x0 +
                           (1.0 / 2.0) * mu *
                               (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) +
                                pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2) - 3));
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
                T x7 = log(x6);
                T x8 = x0 - x3;
                T x9 = mu * x6;
                T x10 = lambda * x7;
                T x11 = -x1 + x4;
                T x12 = x2 - x5;
                T x13 = dx / x6;
                T x14 = f[1] * f[8] - f[2] * f[7];
                T x15 = f[0] * f[8] - f[2] * f[6];
                T x16 = f[0] * f[7] - f[1] * f[6];
                T x17 = f[1] * f[5] - f[2] * f[4];
                T x18 = f[0] * f[5] - f[2] * f[3];
                T x19 = f[0] * f[4] - f[1] * f[3];
                T x20 = mu - x10;
                T x21 = lambda + x20;
                T x22 = x21 * x8;
                T x23 = grad_test[2] * x12;
                T x24 = x11 * x22;
                T x25 = pow(x8, 2);
                T x26 = lambda * x25;
                T x27 = pow(x6, 2);
                T x28 = mu * x27;
                T x29 = x11 * x21;
                T x30 = pow(x12, 2);
                T x31 = pow(x11, 2);
                T x32 = dx / x27;
                T x33 = x14 * x21;
                T x34 = grad_test[0] * x8;
                T x35 = -x33 * x34;
                T x36 = x20 * x6;
                T x37 = f[8] * x36 + x11 * x33;
                T x38 = f[7] * x36 + x12 * x33;
                T x39 = x15 * x21;
                T x40 = grad_test[1] * x11;
                T x41 = -x39 * x40;
                T x42 = f[6] * x36 + x12 * x39;
                T x43 = x6 * (-mu + x10);
                T x44 = f[8] * x43 + x39 * x8;
                T x45 = x16 * x21;
                T x46 = -x23 * x45;
                T x47 = f[6] * x43 + x11 * x45;
                T x48 = f[7] * x43 + x45 * x8;
                T x49 = x17 * x21;
                T x50 = x34 * x49;
                T x51 = f[4] * x36 + x12 * x49;
                T x52 = f[5] * x36 + x11 * x49;
                T x53 = x18 * x21;
                T x54 = x40 * x53;
                T x55 = f[3] * x36 + x12 * x53;
                T x56 = f[5] * x43 + x53 * x8;
                T x57 = x19 * x21;
                T x58 = x23 * x57;
                T x59 = f[4] * x43 + x57 * x8;
                T x60 = f[3] * x43 + x11 * x57;
                T x61 = grad_test[2] * x16;
                T x62 = x15 * x33;
                T x63 = pow(x14, 2);
                T x64 = grad_test[0] * x33;
                T x65 = grad_test[1] * x39;
                T x66 = pow(x16, 2);
                T x67 = pow(x15, 2);
                T x68 = -x17 * x64;
                T x69 = f[2] * x36 + x17 * x39;
                T x70 = f[1] * x36 + x17 * x45;
                T x71 = -x18 * x65;
                T x72 = f[0] * x36 + x18 * x45;
                T x73 = f[2] * x43 + x18 * x33;
                T x74 = grad_test[2] * x19;
                T x75 = -x45 * x74;
                T x76 = f[0] * x43 + x19 * x39;
                T x77 = f[1] * x43 + x19 * x33;
                T x78 = x18 * x49;
                T x79 = pow(x17, 2);
                T x80 = pow(x19, 2);
                T x81 = pow(x18, 2);
                e += dx * ((1.0 / 2.0) * lambda * pow(x7, 2) - mu * x7 +
                           (1.0 / 2.0) * mu *
                               (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) +
                                pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2) - 3));
                lf[0] += x13 * (grad_test[0] * (f[0] * x9 - mu * x8 + x10 * x8) +
                                grad_test[1] * (f[1] * x9 + mu * x11 - x10 * x11) +
                                grad_test[2] * (f[2] * x9 - mu * x12 + x10 * x12));
                lf[1] += x13 * (grad_test[0] * (f[3] * x9 + mu * x14 - x10 * x14) +
                                grad_test[1] * (f[4] * x9 - mu * x15 + x10 * x15) +
                                grad_test[2] * (f[5] * x9 + mu * x16 - x10 * x16));
                lf[2] += x13 * (grad_test[0] * (f[6] * x9 - mu * x17 + x10 * x17) +
                                grad_test[1] * (f[7] * x9 + mu * x18 - x10 * x18) +
                                grad_test[2] * (f[8] * x9 - mu * x19 + x10 * x19));
                bf[0] +=
                    x32 * (grad_trial[0] *
                               (grad_test[0] * (mu * x25 - x26 * x7 + x26 + x28) - grad_test[1] * x24 + x22 * x23) +
                           grad_trial[1] * (-grad_test[0] * x24 +
                                            grad_test[1] * (lambda * x31 + mu * x31 - x10 * x31 + x28) - x23 * x29) +
                           grad_trial[2] * (grad_test[0] * x12 * x22 - grad_test[1] * x12 * x29 +
                                            grad_test[2] * (lambda * x30 + mu * x30 - x10 * x30 + x28)));
                bf[1] += x32 * (grad_trial[0] * (grad_test[1] * x37 - grad_test[2] * x38 + x35) +
                                grad_trial[1] * (grad_test[0] * x44 + grad_test[2] * x42 + x41) +
                                grad_trial[2] * (-grad_test[0] * x48 + grad_test[1] * x47 + x46));
                bf[2] += x32 * (grad_trial[0] * (-grad_test[1] * x52 + grad_test[2] * x51 + x50) +
                                grad_trial[1] * (-grad_test[0] * x56 - grad_test[2] * x55 + x54) +
                                grad_trial[2] * (grad_test[0] * x59 - grad_test[1] * x60 + x58));
                bf[3] += x32 * (grad_trial[0] * (grad_test[1] * x44 - grad_test[2] * x48 + x35) +
                                grad_trial[1] * (grad_test[0] * x37 + grad_test[2] * x47 + x41) +
                                grad_trial[2] * (-grad_test[0] * x38 + grad_test[1] * x42 + x46));
                bf[4] +=
                    x32 * (grad_trial[0] * (grad_test[0] * (lambda * x63 + mu * x63 - x10 * x63 + x28) -
                                            grad_test[1] * x62 + x33 * x61) +
                           grad_trial[1] * (-grad_test[0] * x62 +
                                            grad_test[1] * (lambda * x67 + mu * x67 - x10 * x67 + x28) - x39 * x61) +
                           grad_trial[2] *
                               (grad_test[2] * (lambda * x66 + mu * x66 - x10 * x66 + x28) + x16 * x64 - x16 * x65));
                bf[5] += x32 * (grad_trial[0] * (grad_test[1] * x69 - grad_test[2] * x70 + x68) +
                                grad_trial[1] * (grad_test[0] * x73 + grad_test[2] * x72 + x71) +
                                grad_trial[2] * (-grad_test[0] * x77 + grad_test[1] * x76 + x75));
                bf[6] += x32 * (grad_trial[0] * (-grad_test[1] * x56 + grad_test[2] * x59 + x50) +
                                grad_trial[1] * (-grad_test[0] * x52 - grad_test[2] * x60 + x54) +
                                grad_trial[2] * (grad_test[0] * x51 - grad_test[1] * x55 + x58));
                bf[7] += x32 * (grad_trial[0] * (grad_test[1] * x73 - grad_test[2] * x77 + x68) +
                                grad_trial[1] * (grad_test[0] * x69 + grad_test[2] * x76 + x71) +
                                grad_trial[2] * (-grad_test[0] * x70 + grad_test[1] * x72 + x75));
                bf[8] +=
                    x32 * (grad_trial[0] * (grad_test[0] * (lambda * x79 + mu * x79 - x10 * x79 + x28) -
                                            grad_test[1] * x78 + x49 * x74) +
                           grad_trial[1] * (-grad_test[0] * x78 +
                                            grad_test[1] * (lambda * x81 + mu * x81 - x10 * x81 + x28) - x53 * x74) +
                           grad_trial[2] * (grad_test[0] * x19 * x49 - grad_test[1] * x19 * x53 +
                                            grad_test[2] * (lambda * x80 + mu * x80 - x10 * x80 + x28)));
            }

            T mu{1.0};
            T lambda{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_NeoHookeanOgden_3_IMPL_hpp
