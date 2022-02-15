#ifndef UTOPIA_TPL_HYPERELASTICITY_SaintVenantKirchoff_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_SaintVenantKirchoff_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_SaintVenantKirchoff.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of SaintVenantKirchoff for dimension 2
         */
        template <typename T>
        class SaintVenantKirchoff<T, 2> {
        public:
            static constexpr int Dim = 2;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "SaintVenantKirchoff_2"; }

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    in.get("mu", mu);
                    in.get("lambda", lambda);
                }

                T mu{1.0};
                T lambda{1.0};
            };

            SaintVenantKirchoff(const Params &params) {
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
                T x0 = f[0] * f[1];
                T x1 = f[2] * f[3];
                T x2 = 2 * lambda * x0 + 2 * mu * (2 * x0 + x1);
                T x3 = pow(f[0], 2);
                T x4 = 2 * lambda;
                T x5 = pow(f[1], 2);
                T x6 = pow(f[2], 2);
                T x7 = x5 + x6;
                T x8 = pow(f[3], 2);
                T x9 = x3 + x8;
                T x10 = lambda * (x7 + x9 - 2);
                T x11 = x7 - 1;
                T x12 = 2 * mu;
                T x13 = x9 - 1;
                T x14 = (1.0 / 2.0) * dx;
                T x15 = f[0] * f[3];
                T x16 = f[1] * f[2];
                T x17 = lambda * x16 + mu * x15;
                T x18 = f[0] * f[2];
                T x19 = f[1] * f[3];
                T x20 = grad_test[0] * (lambda * x18 + mu * (2 * x18 + x19));
                T x21 = lambda * x15 + mu * x16;
                T x22 = grad_test[1] * (lambda * x19 + mu * (x18 + 2 * x19));
                T x23 = 2 * lambda * x1 + 2 * mu * (x0 + 2 * x1);
                bf[0] +=
                    x14 * (grad_trial[0] * (grad_test[0] * (x10 + x12 * (x11 + 3 * x3) + x3 * x4) + grad_test[1] * x2) +
                           grad_trial[1] * (grad_test[0] * x2 + grad_test[1] * (x10 + x12 * (x13 + 3 * x5) + x4 * x5)));
                bf[1] += dx * (grad_trial[0] * (grad_test[1] * x17 + x20) + grad_trial[1] * (grad_test[0] * x21 + x22));
                bf[2] += dx * (grad_trial[0] * (grad_test[1] * x21 + x20) + grad_trial[1] * (grad_test[0] * x17 + x22));
                bf[3] += x14 *
                         (grad_trial[0] * (grad_test[0] * (x10 + x12 * (x13 + 3 * x6) + x4 * x6) + grad_test[1] * x23) +
                          grad_trial[1] * (grad_test[0] * x23 + grad_test[1] * (x10 + x12 * (x11 + 3 * x8) + x4 * x8)));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = pow(f[0], 2) + pow(f[2], 2);
                T x1 = pow(f[1], 2) + pow(f[3], 2);
                T x2 = lambda * (x0 + x1 - 2);
                T x3 = f[0] * f[1] + f[2] * f[3];
                T x4 = x0 - 1;
                T x5 = 2 * mu;
                T x6 = x1 - 1;
                T x7 = (1.0 / 2.0) * dx;
                lf[0] += x7 * (grad_test[0] * (f[0] * x2 + x5 * (f[0] * x4 + f[1] * x3)) +
                               grad_test[1] * (f[1] * x2 + x5 * (f[0] * x3 + f[1] * x6)));
                lf[1] += x7 * (grad_test[0] * (f[2] * x2 + x5 * (f[2] * x4 + f[3] * x3)) +
                               grad_test[1] * (f[3] * x2 + x5 * (f[2] * x3 + f[3] * x6)));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = (1.0 / 2.0) * pow(f[0], 2) + (1.0 / 2.0) * pow(f[2], 2);
                T x1 = (1.0 / 2.0) * pow(f[1], 2) + (1.0 / 2.0) * pow(f[3], 2);
                e += dx * ((1.0 / 2.0) * lambda * pow(x0 + x1 - 1, 2) +
                           mu * (pow(x0 - 1.0 / 2.0, 2) + pow(x1 - 1.0 / 2.0, 2) +
                                 2 * pow((1.0 / 2.0) * f[0] * f[1] + (1.0 / 2.0) * f[2] * f[3], 2)));
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
                T x1 = pow(f[2], 2);
                T x2 = (1.0 / 2.0) * x0 + (1.0 / 2.0) * x1;
                T x3 = pow(f[1], 2);
                T x4 = pow(f[3], 2);
                T x5 = (1.0 / 2.0) * x3 + (1.0 / 2.0) * x4;
                T x6 = f[0] * f[1];
                T x7 = f[2] * f[3];
                T x8 = lambda * (x0 + x1 + x3 + x4 - 2);
                T x9 = x6 + x7;
                T x10 = x1 - 1;
                T x11 = x0 + x10;
                T x12 = 2 * mu;
                T x13 = x4 - 1;
                T x14 = x13 + x3;
                T x15 = (1.0 / 2.0) * dx;
                T x16 = 2 * lambda * x6 + 2 * mu * (2 * x6 + x7);
                T x17 = 2 * lambda;
                T x18 = x10 + x3;
                T x19 = x0 + x13;
                T x20 = f[0] * f[3];
                T x21 = f[1] * f[2];
                T x22 = lambda * x21 + mu * x20;
                T x23 = f[0] * f[2];
                T x24 = f[1] * f[3];
                T x25 = grad_test[0] * (lambda * x23 + mu * (2 * x23 + x24));
                T x26 = lambda * x20 + mu * x21;
                T x27 = grad_test[1] * (lambda * x24 + mu * (x23 + 2 * x24));
                T x28 = 2 * lambda * x7 + 2 * mu * (x6 + 2 * x7);
                e += dx * ((1.0 / 2.0) * lambda * pow(x2 + x5 - 1, 2) +
                           mu * (pow(x2 - 1.0 / 2.0, 2) + pow(x5 - 1.0 / 2.0, 2) +
                                 2 * pow((1.0 / 2.0) * x6 + (1.0 / 2.0) * x7, 2)));
                lf[0] += x15 * (grad_test[0] * (f[0] * x8 + x12 * (f[0] * x11 + f[1] * x9)) +
                                grad_test[1] * (f[1] * x8 + x12 * (f[0] * x9 + f[1] * x14)));
                lf[1] += x15 * (grad_test[0] * (f[2] * x8 + x12 * (f[2] * x11 + f[3] * x9)) +
                                grad_test[1] * (f[3] * x8 + x12 * (f[2] * x9 + f[3] * x14)));
                bf[0] += x15 *
                         (grad_trial[0] * (grad_test[0] * (x0 * x17 + x12 * (3 * x0 + x18) + x8) + grad_test[1] * x16) +
                          grad_trial[1] * (grad_test[0] * x16 + grad_test[1] * (x12 * (x19 + 3 * x3) + x17 * x3 + x8)));
                bf[1] += dx * (grad_trial[0] * (grad_test[1] * x22 + x25) + grad_trial[1] * (grad_test[0] * x26 + x27));
                bf[2] += dx * (grad_trial[0] * (grad_test[1] * x26 + x25) + grad_trial[1] * (grad_test[0] * x22 + x27));
                bf[3] += x15 *
                         (grad_trial[0] * (grad_test[0] * (x1 * x17 + x12 * (3 * x1 + x19) + x8) + grad_test[1] * x28) +
                          grad_trial[1] * (grad_test[0] * x28 + grad_test[1] * (x12 * (x18 + 3 * x4) + x17 * x4 + x8)));
            }

            T mu{1.0};
            T lambda{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_SaintVenantKirchoff_2_IMPL_hpp
