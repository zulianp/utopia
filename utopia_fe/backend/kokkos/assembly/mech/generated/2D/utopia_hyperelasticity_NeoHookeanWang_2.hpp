#ifndef UTOPIA_TPL_HYPERELASTICITY_NeoHookeanWang_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeoHookeanWang_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_NeoHookeanWang.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of NeoHookeanWang for dimension 2
         */
        template <typename T>
        class NeoHookeanWang<T, 2> {
        public:
            static constexpr int Dim = 2;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "NeoHookeanWang_2"; }

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
                T x0 = f[2] * f[3];
                T x1 = pow(f[1], 2);
                T x2 = pow(f[2], 2);
                T x3 = x1 + x2;
                T x4 = pow(f[0], 2);
                T x5 = pow(f[3], 2);
                T x6 = x4 + x5;
                T x7 = 5 * x3 + 5 * x6;
                T x8 = f[1] * f[3];
                T x9 = f[0] * f[2];
                T x10 = f[0] * f[3];
                T x11 = f[1] * f[2];
                T x12 = x10 - x11;
                T x13 = 6 * x12;
                T x14 = x0 * x7 + x13 * (x8 - x9);
                T x15 = 9 * pow(x12, 2);
                T x16 = 12 * x12;
                T x17 = -x10 * x16 + x15;
                T x18 = x11 * x16 + x15;
                T x19 = pow(x12, 8.0 / 3.0);
                T x20 = dx / x19;
                T x21 = (1.0 / 9.0) * mu * x20;
                T x22 = f[0] * f[1];
                T x23 = x13 * (-x0 + x22);
                T x24 = 2 * mu;
                T x25 = grad_test[1] * x24 * (x23 + x7 * x9);
                T x26 = 9 * lambda * x19;
                T x27 = x10 * x7;
                T x28 = 3 * x12;
                T x29 = x28 * (x3 + 3 * x4 + 3 * x5);
                T x30 = x7 * x8;
                T x31 = grad_test[0] * x24;
                T x32 = -x24 * (x11 * x7 + x28 * (3 * x1 + 3 * x2 + x6)) + x26;
                T x33 = (1.0 / 18.0) * x20;
                T x34 = x13 * (-x8 + x9) + x22 * x7;
                bf[0] += x21 * (grad_trial[0] * (grad_test[0] * (x17 + x5 * x7) - grad_test[1] * x14) -
                                grad_trial[1] * (grad_test[0] * x14 - grad_test[1] * (x18 + x2 * x7)));
                bf[1] += x33 * (-grad_trial[0] * (grad_test[1] * x32 + x31 * (x13 * (x0 - x22) + x30)) +
                                grad_trial[1] * (grad_test[0] * (-x24 * (-x27 + x29) + x26) - x25));
                bf[2] += -x33 * (grad_trial[0] * (-grad_test[1] * (x24 * (x27 - x29) + x26) + x31 * (-x23 + x30)) +
                                 grad_trial[1] * (grad_test[0] * x32 + x25));
                bf[3] += x21 * (grad_trial[0] * (grad_test[0] * (x1 * x7 + x18) - grad_test[1] * x34) -
                                grad_trial[1] * (grad_test[0] * x34 - grad_test[1] * (x17 + x4 * x7)));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[0] * f[3] - f[1] * f[2];
                T x1 = pow(x0, 5.0 / 3.0);
                T x2 = 3 * lambda * x1;
                T x3 = 3 * x0;
                T x4 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2);
                T x5 = 2 * mu;
                T x6 = (1.0 / 6.0) * dx / x1;
                lf[0] += x6 * (grad_test[0] * (f[3] * x2 + x5 * (f[0] * x3 - f[3] * x4)) -
                               grad_test[1] * (f[2] * x2 - x5 * (f[1] * x3 + f[2] * x4)));
                lf[1] += -x6 * (grad_test[0] * (f[1] * x2 - x5 * (f[1] * x4 + f[2] * x3)) -
                                grad_test[1] * (f[0] * x2 - x5 * (f[0] * x4 - f[3] * x3)));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[0] * f[3] - f[1] * f[2];
                e += dx * ((1.0 / 2.0) * lambda * (x0 - 1) +
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
                T x2 = x0 - x1;
                T x3 = pow(f[1], 2);
                T x4 = pow(f[2], 2);
                T x5 = x3 + x4;
                T x6 = pow(f[0], 2);
                T x7 = pow(f[3], 2);
                T x8 = x6 + x7;
                T x9 = x5 + x8;
                T x10 = pow(x2, 5.0 / 3.0);
                T x11 = 3 * lambda * x10;
                T x12 = 3 * x2;
                T x13 = f[3] * x9;
                T x14 = 2 * mu;
                T x15 = f[2] * x9;
                T x16 = (1.0 / 6.0) * dx / x10;
                T x17 = f[1] * x9;
                T x18 = 5 * x13;
                T x19 = f[1] * f[3];
                T x20 = f[0] * f[2];
                T x21 = 6 * x2;
                T x22 = f[2] * x18 + x21 * (x19 - x20);
                T x23 = 5 * x9;
                T x24 = 9 * pow(x2, 2);
                T x25 = 12 * x2;
                T x26 = -x0 * x25 + x24;
                T x27 = x1 * x25 + x24;
                T x28 = pow(x2, 8.0 / 3.0);
                T x29 = dx / x28;
                T x30 = (1.0 / 9.0) * mu * x29;
                T x31 = 5 * f[0];
                T x32 = f[0] * f[1];
                T x33 = f[2] * f[3];
                T x34 = x21 * (x32 - x33);
                T x35 = grad_test[1] * x14 * (x15 * x31 + x34);
                T x36 = 9 * lambda * x28;
                T x37 = x0 * x23;
                T x38 = x12 * (x5 + 3 * x6 + 3 * x7);
                T x39 = f[1] * x18;
                T x40 = grad_test[0] * x14;
                T x41 = -x14 * (x1 * x23 + x12 * (3 * x3 + 3 * x4 + x8)) + x36;
                T x42 = (1.0 / 18.0) * x29;
                T x43 = x17 * x31 + x21 * (-x19 + x20);
                e += dx * ((1.0 / 2.0) * lambda * (x2 - 1) + (1.0 / 2.0) * mu * (-2 + x9 / pow(x2, 2.0 / 3.0)));
                lf[0] += x16 * (grad_test[0] * (f[3] * x11 + x14 * (f[0] * x12 - x13)) -
                                grad_test[1] * (f[2] * x11 - x14 * (f[1] * x12 + x15)));
                lf[1] += -x16 * (grad_test[0] * (f[1] * x11 - x14 * (f[2] * x12 + x17)) -
                                 grad_test[1] * (f[0] * x11 - x14 * (f[0] * x9 - f[3] * x12)));
                bf[0] += x30 * (grad_trial[0] * (grad_test[0] * (x23 * x7 + x26) - grad_test[1] * x22) -
                                grad_trial[1] * (grad_test[0] * x22 - grad_test[1] * (x23 * x4 + x27)));
                bf[1] += x42 * (-grad_trial[0] * (grad_test[1] * x41 + x40 * (x21 * (-x32 + x33) + x39)) +
                                grad_trial[1] * (grad_test[0] * (-x14 * (-x37 + x38) + x36) - x35));
                bf[2] += -x42 * (grad_trial[0] * (-grad_test[1] * (x14 * (x37 - x38) + x36) + x40 * (-x34 + x39)) +
                                 grad_trial[1] * (grad_test[0] * x41 + x35));
                bf[3] += x30 * (grad_trial[0] * (grad_test[0] * (x23 * x3 + x27) - grad_test[1] * x43) -
                                grad_trial[1] * (grad_test[0] * x43 - grad_test[1] * (x23 * x6 + x26)));
            }

            UTOPIA_FUNCTION void apply(const T *UTOPIA_RESTRICT f,
                                       const T *grad_test,
                                       const T *disp_grad,
                                       const T dx,
                                       T *UTOPIA_RESTRICT res) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[2] * f[3];
                T x1 = pow(f[1], 2);
                T x2 = pow(f[2], 2);
                T x3 = x1 + x2;
                T x4 = pow(f[0], 2);
                T x5 = pow(f[3], 2);
                T x6 = x4 + x5;
                T x7 = 5 * x3 + 5 * x6;
                T x8 = f[1] * f[3];
                T x9 = f[0] * f[2];
                T x10 = f[0] * f[3];
                T x11 = f[1] * f[2];
                T x12 = x10 - x11;
                T x13 = 6 * x12;
                T x14 = x0 * x7 + x13 * (x8 - x9);
                T x15 = 2 * mu;
                T x16 = disp_grad[1] * x15;
                T x17 = f[0] * f[1];
                T x18 = x13 * (x0 - x17) + x7 * x8;
                T x19 = disp_grad[2] * x15;
                T x20 = 9 * pow(x12, 2);
                T x21 = 12 * x12;
                T x22 = -x10 * x21 + x20;
                T x23 = disp_grad[0] * x15;
                T x24 = pow(x12, 8.0 / 3.0);
                T x25 = 9 * lambda * x24;
                T x26 = 3 * x12;
                T x27 = -x15 * (-x10 * x7 + x26 * (x3 + 3 * x4 + 3 * x5)) + x25;
                T x28 = x13 * (-x0 + x17) + x7 * x9;
                T x29 = disp_grad[3] * x15;
                T x30 = x11 * x21 + x20;
                T x31 = -x15 * (x11 * x7 + x26 * (3 * x1 + 3 * x2 + x6)) + x25;
                T x32 = (1.0 / 18.0) * dx / x24;
                T x33 = x13 * (-x8 + x9) + x17 * x7;
                res[0] += x32 * (grad_test[0] * (disp_grad[3] * x27 - x14 * x16 - x18 * x19 + x23 * (x22 + x5 * x7)) -
                                 grad_test[1] * (disp_grad[2] * x31 + x14 * x23 - x16 * (x2 * x7 + x30) + x28 * x29));
                res[1] += x32 * (-grad_test[0] * (disp_grad[1] * x31 + x18 * x23 - x19 * (x1 * x7 + x30) + x29 * x33) +
                                 grad_test[1] * (disp_grad[0] * x27 - x16 * x28 - x19 * x33 + x29 * (x22 + x4 * x7)));
            }

            T mu{1.0};
            T lambda{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_NeoHookeanWang_2_IMPL_hpp
