#ifndef UTOPIA_TPL_HYPERELASTICITY_NeoHookeanBower_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeoHookeanBower_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_NeoHookeanBower.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of NeoHookeanBower for dimension 2
         */
        template <typename T>
        class NeoHookeanBower<T, 2> {
        public:
            static constexpr int Dim = 2;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "NeoHookeanBower_2"; }

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

            NeoHookeanBower(const Params &params) {
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
                T x0 = f[2] * f[3];
                T x1 = f[0] * f[3];
                T x2 = f[1] * f[2];
                T x3 = -x2;
                T x4 = x1 + x3;
                T x5 = pow(x4, 8.0 / 3.0);
                T x6 = 9 * lambda;
                T x7 = x5 * x6;
                T x8 = x0 * x7;
                T x9 = pow(f[0], 2);
                T x10 = pow(f[3], 2);
                T x11 = x10 + x9;
                T x12 = pow(f[1], 2);
                T x13 = pow(f[2], 2);
                T x14 = x12 + x13;
                T x15 = 5 * x11 + 5 * x14;
                T x16 = x0 * x15;
                T x17 = f[1] * f[3];
                T x18 = f[0] * f[2];
                T x19 = 6 * x4;
                T x20 = x19 * (x17 - x18);
                T x21 = pow(x4, 7.0 / 3.0);
                T x22 = grad_test[1] * x21;
                T x23 = pow(x4, 5);
                T x24 = x23 * x6;
                T x25 = 9 * pow(x4, 2);
                T x26 = 12 * x4;
                T x27 = -x1 * x26 + x25;
                T x28 = mu * x21;
                T x29 = x19 * (-x17 + x18);
                T x30 = grad_test[0] * x21;
                T x31 = x2 * x26 + x25;
                T x32 = (1.0 / 9.0) * dx;
                T x33 = x32 / x23;
                T x34 = f[0] * f[1];
                T x35 = x19 * (-x0 + x34);
                T x36 = grad_test[1] * (mu * (x15 * x18 + x35) + x18 * x7);
                T x37 = 2 * x1;
                T x38 = 3 * x4;
                T x39 = mu * (-x1 * x15 + x38 * (3 * x10 + x14 + 3 * x9));
                T x40 = x17 * x7;
                T x41 = x15 * x17;
                T x42 = mu * (x15 * x2 + x38 * (x11 + 3 * x12 + 3 * x13)) + x7 * (-x1 + 2 * x2 + 1);
                T x43 = x32 / x5;
                T x44 = x34 * x7;
                T x45 = x15 * x34;
                bf[0] += x33 * (grad_trial[0] * (grad_test[0] * (x10 * x24 + x28 * (x10 * x15 + x27)) -
                                                 x22 * (mu * (x16 + x20) + x8)) -
                                grad_trial[1] * (-grad_test[1] * (x13 * x24 + x28 * (x13 * x15 + x31)) +
                                                 x30 * (-mu * (-x16 + x29) + x8)));
                bf[1] += -x43 * (grad_trial[0] * (grad_test[0] * (mu * (-x35 + x41) + x40) - grad_test[1] * x42) +
                                 grad_trial[1] * (grad_test[0] * (x39 + x7 * (x2 - x37 + 1)) + x36));
                bf[2] += x43 * (-grad_trial[0] * (grad_test[0] * (-mu * (x35 - x41) + x40) -
                                                  grad_test[1] * (-x39 + x7 * (x3 + x37 - 1))) +
                                grad_trial[1] * (grad_test[0] * x42 - x36));
                bf[3] +=
                    x33 * (grad_trial[0] *
                               (grad_test[0] * (x12 * x24 + x28 * (x12 * x15 + x31)) - x22 * (mu * (x29 + x45) + x44)) -
                           grad_trial[1] *
                               (-grad_test[1] * (x24 * x9 + x28 * (x15 * x9 + x27)) + x30 * (-mu * (x20 - x45) + x44)));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[0] * f[3];
                T x1 = f[1] * f[2];
                T x2 = x0 - x1;
                T x3 = pow(x2, 5.0 / 3.0);
                T x4 = 3 * lambda * x3 * (-x0 + x1 + 1);
                T x5 = 3 * x2;
                T x6 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2);
                T x7 = (1.0 / 3.0) * dx / x3;
                lf[0] += x7 * (-grad_test[0] * (f[3] * x4 - mu * (f[0] * x5 - f[3] * x6)) +
                               grad_test[1] * (f[2] * x4 + mu * (f[1] * x5 + f[2] * x6)));
                lf[1] += -x7 * (-grad_test[0] * (f[1] * x4 + mu * (f[1] * x6 + f[2] * x5)) +
                                grad_test[1] * (f[0] * x4 + mu * (f[0] * x6 - f[3] * x5)));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[0] * f[3] - f[1] * f[2];
                e += dx * ((1.0 / 2.0) * lambda * pow(x0 - 1, 2) +
                           (1.0 / 2.0) * mu *
                               (-2 + (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2)) / pow(x0, 2.0 / 3.0)));
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
                T x0 = f[0] * f[3];
                T x1 = f[1] * f[2];
                T x2 = -x1;
                T x3 = x0 + x2;
                T x4 = pow(f[0], 2);
                T x5 = pow(f[3], 2);
                T x6 = x4 + x5;
                T x7 = pow(f[1], 2);
                T x8 = pow(f[2], 2);
                T x9 = x7 + x8;
                T x10 = x6 + x9;
                T x11 = pow(x3, 5.0 / 3.0);
                T x12 = 1 - x0;
                T x13 = 3 * lambda * x11 * (x1 + x12);
                T x14 = 3 * x3;
                T x15 = f[2] * x10;
                T x16 = f[3] * x10;
                T x17 = (1.0 / 3.0) * dx / x11;
                T x18 = f[0] * x10;
                T x19 = f[2] * f[3];
                T x20 = pow(x3, 8.0 / 3.0);
                T x21 = 9 * lambda;
                T x22 = x20 * x21;
                T x23 = x19 * x22;
                T x24 = 5 * x15;
                T x25 = f[3] * x24;
                T x26 = f[1] * f[3];
                T x27 = f[0] * f[2];
                T x28 = 6 * x3;
                T x29 = x28 * (x26 - x27);
                T x30 = pow(x3, 7.0 / 3.0);
                T x31 = grad_test[1] * x30;
                T x32 = pow(x3, 5);
                T x33 = x21 * x32;
                T x34 = 5 * x10;
                T x35 = 9 * pow(x3, 2);
                T x36 = 12 * x3;
                T x37 = -x0 * x36 + x35;
                T x38 = mu * x30;
                T x39 = x28 * (-x26 + x27);
                T x40 = grad_test[0] * x30;
                T x41 = x1 * x36 + x35;
                T x42 = (1.0 / 9.0) * dx;
                T x43 = x42 / x32;
                T x44 = f[0] * f[1];
                T x45 = x28 * (-x19 + x44);
                T x46 = grad_test[1] * (mu * (f[0] * x24 + x45) + x22 * x27);
                T x47 = 2 * x0;
                T x48 = mu * (-x0 * x34 + x14 * (3 * x4 + 3 * x5 + x9));
                T x49 = x22 * x26;
                T x50 = 5 * f[1];
                T x51 = x16 * x50;
                T x52 = mu * (x1 * x34 + x14 * (x6 + 3 * x7 + 3 * x8)) + x22 * (2 * x1 + x12);
                T x53 = x42 / x20;
                T x54 = x22 * x44;
                T x55 = x18 * x50;
                e += dx * ((1.0 / 2.0) * lambda * pow(x3 - 1, 2) + (1.0 / 2.0) * mu * (x10 / pow(x3, 2.0 / 3.0) - 2));
                lf[0] += x17 * (-grad_test[0] * (f[3] * x13 - mu * (f[0] * x14 - x16)) +
                                grad_test[1] * (f[2] * x13 + mu * (f[1] * x14 + x15)));
                lf[1] += -x17 * (-grad_test[0] * (f[1] * x13 + mu * (f[1] * x10 + f[2] * x14)) +
                                 grad_test[1] * (f[0] * x13 + mu * (-f[3] * x14 + x18)));
                bf[0] += x43 * (grad_trial[0] * (grad_test[0] * (x33 * x5 + x38 * (x34 * x5 + x37)) -
                                                 x31 * (mu * (x25 + x29) + x23)) -
                                grad_trial[1] * (-grad_test[1] * (x33 * x8 + x38 * (x34 * x8 + x41)) +
                                                 x40 * (-mu * (-x25 + x39) + x23)));
                bf[1] += -x53 * (grad_trial[0] * (grad_test[0] * (mu * (-x45 + x51) + x49) - grad_test[1] * x52) +
                                 grad_trial[1] * (grad_test[0] * (x22 * (x1 - x47 + 1) + x48) + x46));
                bf[2] += x53 * (-grad_trial[0] * (grad_test[0] * (-mu * (x45 - x51) + x49) -
                                                  grad_test[1] * (x22 * (x2 + x47 - 1) - x48)) +
                                grad_trial[1] * (grad_test[0] * x52 - x46));
                bf[3] +=
                    x43 * (grad_trial[0] *
                               (grad_test[0] * (x33 * x7 + x38 * (x34 * x7 + x41)) - x31 * (mu * (x39 + x55) + x54)) -
                           grad_trial[1] *
                               (-grad_test[1] * (x33 * x4 + x38 * (x34 * x4 + x37)) + x40 * (-mu * (x29 - x55) + x54)));
            }

            T mu{1.0};
            T lambda{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_NeoHookeanBower_2_IMPL_hpp
