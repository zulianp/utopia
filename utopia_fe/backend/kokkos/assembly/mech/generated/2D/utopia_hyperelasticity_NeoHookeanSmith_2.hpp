#ifndef UTOPIA_TPL_HYPERELASTICITY_NeoHookeanSmith_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeoHookeanSmith_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_NeoHookeanSmith.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of NeoHookeanSmith for dimension 2
         */
        template <typename T>
        class NeoHookeanSmith<T, 2> {
        public:
            static constexpr int Dim = 2;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "NeoHookeanSmith_2"; }

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
                T x0 = 16 * mu;
                T x1 = f[0] * x0;
                T x2 = pow(f[0], 2);
                T x3 = pow(f[1], 2);
                T x4 = pow(f[2], 2);
                T x5 = pow(f[3], 2);
                T x6 = x2 + x3 + x4 + x5 + 1;
                T x7 = pow(x6, 2);
                T x8 = 11 * lambda;
                T x9 = x7 * x8;
                T x10 = f[3] * x9;
                T x11 = f[1] * x1 - f[2] * x10;
                T x12 = 8 * mu;
                T x13 = -x12 * x6;
                T x14 = (1.0 / 6.0) * dx / x7;
                T x15 = grad_test[0] * (-f[1] * x10 + f[2] * x1);
                T x16 = f[1] * f[2];
                T x17 = 6 * mu;
                T x18 = f[0] * f[3];
                T x19 = x8 * (x16 - x18 + 1);
                T x20 = x0 * x16 + x7 * (x16 * x8 + x17 + x19);
                T x21 = f[3] * x0;
                T x22 = f[0] * x9;
                T x23 = -grad_test[1] * (-f[1] * x21 + f[2] * x22);
                T x24 = x0 * x18 + x7 * (-x17 + x18 * x8 - x19);
                T x25 = f[1] * x22 - f[2] * x21;
                bf[0] += x14 *
                         (grad_trial[0] * (grad_test[0] * (x0 * x2 + x13 + x7 * (x12 + x5 * x8)) + grad_test[1] * x11) +
                          grad_trial[1] * (grad_test[0] * x11 + grad_test[1] * (x0 * x3 + x13 + x7 * (x12 + x4 * x8))));
                bf[1] +=
                    x14 * (grad_trial[0] * (grad_test[1] * x20 + x15) + grad_trial[1] * (grad_test[0] * x24 + x23));
                bf[2] +=
                    x14 * (grad_trial[0] * (grad_test[1] * x24 + x15) + grad_trial[1] * (grad_test[0] * x20 + x23));
                bf[3] += x14 *
                         (grad_trial[0] * (grad_test[0] * (x0 * x4 + x13 + x7 * (x12 + x3 * x8)) - grad_test[1] * x25) -
                          grad_trial[1] * (grad_test[0] * x25 - grad_test[1] * (x0 * x5 + x13 + x7 * (x12 + x2 * x8))));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = 8 * mu;
                T x1 = f[1] * x0;
                T x2 = pow(f[0], 2);
                T x3 = pow(f[1], 2);
                T x4 = pow(f[2], 2);
                T x5 = pow(f[3], 2);
                T x6 = x2 + x3 + x4 + x5 + 1;
                T x7 = 6 * mu;
                T x8 = 11 * lambda * (-f[0] * f[3] + f[1] * f[2] + 1);
                T x9 = x7 + x8;
                T x10 = f[0] * x0;
                T x11 = dx / (6 * x2 + 6 * x3 + 6 * x4 + 6 * x5 + 6);
                T x12 = f[3] * x0;
                T x13 = f[2] * x0;
                lf[0] += -x11 *
                         (grad_test[0] * (x10 - x6 * (-f[3] * x9 + x10)) + grad_test[1] * (x1 - x6 * (f[2] * x9 + x1)));
                lf[1] += -x11 * (grad_test[0] * (x13 - x6 * (-f[1] * (-x7 - x8) + x13)) +
                                 grad_test[1] * (x12 + x6 * (f[0] * x9 - x12)));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2);
                T x1 = (2.0 / 3.0) * mu;
                e += dx * ((11.0 / 12.0) * lambda * pow(f[0] * f[3] - f[1] * f[2] - 1 - 6.0 / 11.0 * mu / lambda, 2) +
                           x1 * (x0 - 2) - x1 * log(x0 + 1));
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
                T x4 = x0 + x1 + x2 + x3;
                T x5 = (2.0 / 3.0) * mu;
                T x6 = x4 + 1;
                T x7 = f[0] * f[3];
                T x8 = f[1] * f[2];
                T x9 = 8 * mu;
                T x10 = f[1] * x9;
                T x11 = 6 * mu;
                T x12 = 11 * lambda;
                T x13 = x12 * (-x7 + x8 + 1);
                T x14 = x11 + x13;
                T x15 = f[0] * x9;
                T x16 = dx / (6 * x0 + 6 * x1 + 6 * x2 + 6 * x3 + 6);
                T x17 = f[3] * x9;
                T x18 = f[2] * x9;
                T x19 = -x11 - x13;
                T x20 = 16 * mu;
                T x21 = f[0] * x20;
                T x22 = pow(x6, 2);
                T x23 = x12 * x22;
                T x24 = f[3] * x23;
                T x25 = f[1] * x21 - f[2] * x24;
                T x26 = -x6 * x9;
                T x27 = (1.0 / 6.0) * dx / x22;
                T x28 = grad_test[0] * (-f[1] * x24 + f[2] * x21);
                T x29 = x20 * x8 + x22 * (x12 * x8 + x14);
                T x30 = f[3] * x20;
                T x31 = f[0] * x23;
                T x32 = -grad_test[1] * (-f[1] * x30 + f[2] * x31);
                T x33 = x20 * x7 + x22 * (x12 * x7 + x19);
                T x34 = f[1] * x31 - f[2] * x30;
                e += dx * ((11.0 / 12.0) * lambda * pow(x7 - x8 - 1 - 6.0 / 11.0 * mu / lambda, 2) + x5 * (x4 - 2) -
                           x5 * log(x6));
                lf[0] += -x16 * (grad_test[0] * (x15 - x6 * (-f[3] * x14 + x15)) +
                                 grad_test[1] * (x10 - x6 * (f[2] * x14 + x10)));
                lf[1] += -x16 * (grad_test[0] * (x18 - x6 * (-f[1] * x19 + x18)) +
                                 grad_test[1] * (x17 + x6 * (f[0] * x14 - x17)));
                bf[0] +=
                    x27 *
                    (grad_trial[0] * (grad_test[0] * (x0 * x20 + x22 * (x12 * x3 + x9) + x26) + grad_test[1] * x25) +
                     grad_trial[1] * (grad_test[0] * x25 + grad_test[1] * (x1 * x20 + x22 * (x12 * x2 + x9) + x26)));
                bf[1] +=
                    x27 * (grad_trial[0] * (grad_test[1] * x29 + x28) + grad_trial[1] * (grad_test[0] * x33 + x32));
                bf[2] +=
                    x27 * (grad_trial[0] * (grad_test[1] * x33 + x28) + grad_trial[1] * (grad_test[0] * x29 + x32));
                bf[3] +=
                    x27 *
                    (grad_trial[0] * (grad_test[0] * (x2 * x20 + x22 * (x1 * x12 + x9) + x26) - grad_test[1] * x34) -
                     grad_trial[1] * (grad_test[0] * x34 - grad_test[1] * (x20 * x3 + x22 * (x0 * x12 + x9) + x26)));
            }

            UTOPIA_FUNCTION void apply(const T *UTOPIA_RESTRICT f,
                                       const T *grad_test,
                                       const T *disp_grad,
                                       const T dx,
                                       T *UTOPIA_RESTRICT res) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = 16 * mu;
                T x1 = f[0] * x0;
                T x2 = pow(f[0], 2);
                T x3 = pow(f[1], 2);
                T x4 = pow(f[2], 2);
                T x5 = pow(f[3], 2);
                T x6 = x2 + x3 + x4 + x5 + 1;
                T x7 = pow(x6, 2);
                T x8 = 11 * lambda;
                T x9 = x7 * x8;
                T x10 = f[3] * x9;
                T x11 = f[1] * x1 - f[2] * x10;
                T x12 = -f[1] * x10 + f[2] * x1;
                T x13 = f[0] * f[3];
                T x14 = 6 * mu;
                T x15 = f[1] * f[2];
                T x16 = x8 * (-x13 + x15 + 1);
                T x17 = x0 * x13 + x7 * (x13 * x8 - x14 - x16);
                T x18 = 8 * mu;
                T x19 = -x18 * x6;
                T x20 = f[3] * x0;
                T x21 = f[0] * x9;
                T x22 = -f[1] * x20 + f[2] * x21;
                T x23 = x0 * x15 + x7 * (x14 + x15 * x8 + x16);
                T x24 = (1.0 / 6.0) * dx / x7;
                T x25 = f[1] * x21 - f[2] * x20;
                res[0] +=
                    x24 * (grad_test[0] * (disp_grad[0] * (x0 * x2 + x19 + x7 * (x18 + x5 * x8)) + disp_grad[1] * x11 +
                                           disp_grad[2] * x12 + disp_grad[3] * x17) +
                           grad_test[1] * (disp_grad[0] * x11 + disp_grad[1] * (x0 * x3 + x19 + x7 * (x18 + x4 * x8)) +
                                           disp_grad[2] * x23 - disp_grad[3] * x22));
                res[1] +=
                    x24 * (grad_test[0] * (disp_grad[0] * x12 + disp_grad[1] * x23 +
                                           disp_grad[2] * (x0 * x4 + x19 + x7 * (x18 + x3 * x8)) - disp_grad[3] * x25) +
                           grad_test[1] * (disp_grad[0] * x17 - disp_grad[1] * x22 - disp_grad[2] * x25 +
                                           disp_grad[3] * (x0 * x5 + x19 + x7 * (x18 + x2 * x8))));
            }

            T mu{1.0};
            T lambda{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_NeoHookeanSmith_2_IMPL_hpp
