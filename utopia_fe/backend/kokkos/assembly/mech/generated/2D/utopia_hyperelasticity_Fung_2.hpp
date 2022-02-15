#ifndef UTOPIA_TPL_HYPERELASTICITY_Fung_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_Fung_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_Fung.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Fung for dimension 2
         */
        template <typename T>
        class Fung<T, 2> {
        public:
            static constexpr int Dim = 2;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Fung_2"; }

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    in.get("a", a);
                    in.get("b", b);
                    in.get("c", c);
                    in.get("k", k);
                }

                T a{1.0};
                T b{1.0};
                T c{1.0};
                T k{1};
            };

            Fung(const Params &params) {
                a = params.a;
                b = params.b;
                c = params.c;
                k = params.k;
            }

            UTOPIA_FUNCTION void hessian(const T *UTOPIA_RESTRICT f,
                                         const T *grad_test,
                                         const T *grad_trial,
                                         const T dx,
                                         T *UTOPIA_RESTRICT bf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = 2 * k;
                T x1 = f[3] * x0;
                T x2 = pow(f[0], 2);
                T x3 = pow(f[1], 2);
                T x4 = pow(f[2], 2);
                T x5 = pow(f[3], 2);
                T x6 = c * exp((1.0 / 2.0) * b * (x2 + x3 + x4 + x5 - 2));
                T x7 = pow(b, 2) * x6;
                T x8 = f[0] * x7;
                T x9 = f[1] * x8 - f[2] * x1;
                T x10 = a + b * x6;
                T x11 = (1.0 / 2.0) * dx;
                T x12 = grad_test[0] * (-f[1] * x1 + f[2] * x8);
                T x13 = f[1] * f[2];
                T x14 = f[0] * f[3];
                T x15 = x0 * (x13 - x14 + 1);
                T x16 = x0 * x13 + x13 * x7 + x15;
                T x17 = f[0] * x0;
                T x18 = f[3] * x7;
                T x19 = grad_test[1] * (f[1] * x18 - f[2] * x17);
                T x20 = x0 * x14 + x14 * x7 - x15;
                T x21 = -f[1] * x17 + f[2] * x18;
                bf[0] += x11 * (grad_trial[0] * (grad_test[0] * (x0 * x5 + x10 + x2 * x7) + grad_test[1] * x9) +
                                grad_trial[1] * (grad_test[0] * x9 + grad_test[1] * (x0 * x4 + x10 + x3 * x7)));
                bf[1] +=
                    x11 * (grad_trial[0] * (grad_test[1] * x16 + x12) + grad_trial[1] * (grad_test[0] * x20 + x19));
                bf[2] +=
                    x11 * (grad_trial[0] * (grad_test[1] * x20 + x12) + grad_trial[1] * (grad_test[0] * x16 + x19));
                bf[3] += x11 * (grad_trial[0] * (grad_test[0] * (x0 * x3 + x10 + x4 * x7) + grad_test[1] * x21) +
                                grad_trial[1] * (grad_test[0] * x21 + grad_test[1] * (x0 * x2 + x10 + x5 * x7)));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = 2 * k * (-f[0] * f[3] + f[1] * f[2] + 1);
                T x1 = b * c * exp((1.0 / 2.0) * b * (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) - 2));
                T x2 = (1.0 / 2.0) * dx;
                lf[0] += x2 * (grad_test[0] * (a * f[0] + f[0] * x1 - f[3] * x0) +
                               grad_test[1] * (a * f[1] + f[1] * x1 + f[2] * x0));
                lf[1] += x2 * (grad_test[0] * (a * f[2] + f[1] * x0 + f[2] * x1) +
                               grad_test[1] * (a * f[3] - f[0] * x0 + f[3] * x1));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = (1.0 / 2.0) * pow(f[0], 2) + (1.0 / 2.0) * pow(f[1], 2) + (1.0 / 2.0) * pow(f[2], 2) +
                       (1.0 / 2.0) * pow(f[3], 2) - 1;
                e += dx * ((1.0 / 2.0) * a * x0 + (1.0 / 2.0) * c * (exp(b * x0) - 1) +
                           (1.0 / 2.0) * k * pow(f[0] * f[3] - f[1] * f[2] - 1, 2));
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
                T x2 = pow(f[0], 2);
                T x3 = pow(f[1], 2);
                T x4 = pow(f[2], 2);
                T x5 = pow(f[3], 2);
                T x6 = (1.0 / 2.0) * x2 + (1.0 / 2.0) * x3 + (1.0 / 2.0) * x4 + (1.0 / 2.0) * x5 - 1;
                T x7 = 2 * k;
                T x8 = x7 * (-x0 + x1 + 1);
                T x9 = c * exp((1.0 / 2.0) * b * (x2 + x3 + x4 + x5 - 2));
                T x10 = b * x9;
                T x11 = (1.0 / 2.0) * dx;
                T x12 = f[3] * x7;
                T x13 = pow(b, 2) * x9;
                T x14 = f[0] * x13;
                T x15 = f[1] * x14 - f[2] * x12;
                T x16 = a + x10;
                T x17 = grad_test[0] * (-f[1] * x12 + f[2] * x14);
                T x18 = x1 * x13 + x1 * x7 + x8;
                T x19 = f[0] * x7;
                T x20 = f[3] * x13;
                T x21 = grad_test[1] * (f[1] * x20 - f[2] * x19);
                T x22 = x0 * x13 + x0 * x7 - x8;
                T x23 = -f[1] * x19 + f[2] * x20;
                e += dx * ((1.0 / 2.0) * a * x6 + (1.0 / 2.0) * c * (exp(b * x6) - 1) +
                           (1.0 / 2.0) * k * pow(x0 - x1 - 1, 2));
                lf[0] += x11 * (grad_test[0] * (a * f[0] + f[0] * x10 - f[3] * x8) +
                                grad_test[1] * (a * f[1] + f[1] * x10 + f[2] * x8));
                lf[1] += x11 * (grad_test[0] * (a * f[2] + f[1] * x8 + f[2] * x10) +
                                grad_test[1] * (a * f[3] - f[0] * x8 + f[3] * x10));
                bf[0] += x11 * (grad_trial[0] * (grad_test[0] * (x13 * x2 + x16 + x5 * x7) + grad_test[1] * x15) +
                                grad_trial[1] * (grad_test[0] * x15 + grad_test[1] * (x13 * x3 + x16 + x4 * x7)));
                bf[1] +=
                    x11 * (grad_trial[0] * (grad_test[1] * x18 + x17) + grad_trial[1] * (grad_test[0] * x22 + x21));
                bf[2] +=
                    x11 * (grad_trial[0] * (grad_test[1] * x22 + x17) + grad_trial[1] * (grad_test[0] * x18 + x21));
                bf[3] += x11 * (grad_trial[0] * (grad_test[0] * (x13 * x4 + x16 + x3 * x7) + grad_test[1] * x23) +
                                grad_trial[1] * (grad_test[0] * x23 + grad_test[1] * (x13 * x5 + x16 + x2 * x7)));
            }

            T a{1.0};
            T b{1.0};
            T c{1.0};
            T k{1};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_Fung_2_IMPL_hpp
