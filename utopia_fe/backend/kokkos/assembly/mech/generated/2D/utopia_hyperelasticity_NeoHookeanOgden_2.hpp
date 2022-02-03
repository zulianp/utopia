#ifndef UTOPIA_TPL_HYPERELASTICITY_NeoHookeanOgden_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_NeoHookeanOgden_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_NeoHookeanOgden.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of NeoHookeanOgden for dimension 2
         */
        template <typename T>
        class NeoHookeanOgden<T, 2> {
        public:
            static constexpr int Dim = 2;

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
                T x0 = f[0] * f[3];
                T x1 = f[1] * f[2];
                T x2 = x0 - x1;
                T x3 = log(x2);
                T x4 = lambda * x3;
                T x5 = mu - x4;
                T x6 = lambda + x5;
                T x7 = f[3] * x6;
                T x8 = f[2] * x7;
                T x9 = pow(f[3], 2);
                T x10 = lambda * x9;
                T x11 = pow(x2, 2);
                T x12 = mu * x11;
                T x13 = pow(f[2], 2);
                T x14 = dx / x11;
                T x15 = f[1] * grad_test[0];
                T x16 = x15 * x7;
                T x17 = x1 * x6 + x2 * x5;
                T x18 = f[0] * x6;
                T x19 = grad_test[1] * x18;
                T x20 = -f[2] * x19;
                T x21 = x0 * x6 + x2 * (-mu + x4);
                T x22 = pow(f[1], 2);
                T x23 = pow(f[0], 2);
                bf[0] += x14 * (grad_trial[0] * (grad_test[0] * (mu * x9 - x10 * x3 + x10 + x12) - grad_test[1] * x8) -
                                grad_trial[1] *
                                    (grad_test[0] * x8 - grad_test[1] * (lambda * x13 + mu * x13 + x12 - x13 * x4)));
                bf[1] +=
                    x14 * (-grad_trial[0] * (-grad_test[1] * x17 + x16) + grad_trial[1] * (grad_test[0] * x21 + x20));
                bf[2] +=
                    x14 * (-grad_trial[0] * (-grad_test[1] * x21 + x16) + grad_trial[1] * (grad_test[0] * x17 + x20));
                bf[3] +=
                    x14 * (grad_trial[0] * (-f[1] * x19 + grad_test[0] * (lambda * x22 + mu * x22 + x12 - x22 * x4)) -
                           grad_trial[1] * (-grad_test[1] * (lambda * x23 + mu * x23 + x12 - x23 * x4) + x15 * x18));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[3] * mu;
                T x1 = f[0] * f[3] - f[1] * f[2];
                T x2 = f[0] * mu;
                T x3 = lambda * log(x1);
                T x4 = f[2] * mu;
                T x5 = f[1] * mu;
                T x6 = dx / x1;
                lf[0] += x6 * (grad_test[0] * (f[3] * x3 - x0 + x1 * x2) + grad_test[1] * (-f[2] * x3 + x1 * x5 + x4));
                lf[1] += x6 * (grad_test[0] * (-f[1] * x3 + x1 * x4 + x5) + grad_test[1] * (f[0] * x3 + x0 * x1 - x2));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = log(f[0] * f[3] - f[1] * f[2]);
                e += dx * ((1.0 / 2.0) * lambda * pow(x0, 2) - mu * x0 +
                           (1.0 / 2.0) * mu * (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) - 2));
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
                T x3 = log(x2);
                T x4 = pow(f[0], 2);
                T x5 = pow(f[1], 2);
                T x6 = pow(f[2], 2);
                T x7 = pow(f[3], 2);
                T x8 = f[3] * mu;
                T x9 = f[0] * mu;
                T x10 = lambda * x3;
                T x11 = f[2] * mu;
                T x12 = f[1] * mu;
                T x13 = dx / x2;
                T x14 = mu - x10;
                T x15 = lambda + x14;
                T x16 = f[3] * x15;
                T x17 = f[2] * x16;
                T x18 = lambda * x7;
                T x19 = pow(x2, 2);
                T x20 = mu * x19;
                T x21 = dx / x19;
                T x22 = f[1] * grad_test[0];
                T x23 = x16 * x22;
                T x24 = x1 * x15 + x14 * x2;
                T x25 = f[0] * x15;
                T x26 = grad_test[1] * x25;
                T x27 = -f[2] * x26;
                T x28 = x0 * x15 + x2 * (-mu + x10);
                e += dx * ((1.0 / 2.0) * lambda * pow(x3, 2) - mu * x3 + (1.0 / 2.0) * mu * (x4 + x5 + x6 + x7 - 2));
                lf[0] +=
                    x13 * (grad_test[0] * (f[3] * x10 + x2 * x9 - x8) + grad_test[1] * (-f[2] * x10 + x11 + x12 * x2));
                lf[1] +=
                    x13 * (grad_test[0] * (-f[1] * x10 + x11 * x2 + x12) + grad_test[1] * (f[0] * x10 + x2 * x8 - x9));
                bf[0] += x21 * (grad_trial[0] * (grad_test[0] * (mu * x7 - x18 * x3 + x18 + x20) - grad_test[1] * x17) -
                                grad_trial[1] *
                                    (grad_test[0] * x17 - grad_test[1] * (lambda * x6 + mu * x6 - x10 * x6 + x20)));
                bf[1] +=
                    x21 * (-grad_trial[0] * (-grad_test[1] * x24 + x23) + grad_trial[1] * (grad_test[0] * x28 + x27));
                bf[2] +=
                    x21 * (-grad_trial[0] * (-grad_test[1] * x28 + x23) + grad_trial[1] * (grad_test[0] * x24 + x27));
                bf[3] +=
                    x21 * (grad_trial[0] * (-f[1] * x26 + grad_test[0] * (lambda * x5 + mu * x5 - x10 * x5 + x20)) -
                           grad_trial[1] * (-grad_test[1] * (lambda * x4 + mu * x4 - x10 * x4 + x20) + x22 * x25));
            }

            T mu{1.0};
            T lambda{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_NeoHookeanOgden_2_IMPL_hpp
