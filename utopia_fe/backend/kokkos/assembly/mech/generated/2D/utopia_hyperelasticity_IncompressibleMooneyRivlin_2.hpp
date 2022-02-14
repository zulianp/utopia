#ifndef UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_IncompressibleMooneyRivlin.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of IncompressibleMooneyRivlin for dimension 2
         */
        template <typename T>
        class IncompressibleMooneyRivlin<T, 2> {
        public:
            static constexpr int Dim = 2;

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    in.get("C1", C1);
                    in.get("C2", C2);
                }

                T C1{1.0};
                T C2{1.0};
            };

            IncompressibleMooneyRivlin(const Params &params) {
                C1 = params.C1;
                C2 = params.C2;
            }

            UTOPIA_FUNCTION void hessian(const T *UTOPIA_RESTRICT f,
                                         const T p,
                                         const T *grad_test,
                                         const T *grad_trial,
                                         const T fun_test,
                                         const T fun_trial,
                                         const T dx,
                                         T *UTOPIA_RESTRICT bf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[2] * grad_test[1];
                T x1 = f[3] * grad_test[0];
                T x2 = 2 * dx;
                T x3 = 2 * C2;
                T x4 = f[1] * x1 * x3;
                T x5 = f[0] * f[3];
                T x6 = f[1] * f[2];
                T x7 = p + x3 * (x5 - 2 * x6);
                T x8 = f[0] * x0 * x3;
                T x9 = p + x3 * (2 * x5 - x6);
                T x10 = f[0] * grad_test[1];
                T x11 = f[1] * grad_test[0];
                T x12 = dx * fun_test;
                T x13 = dx * fun_trial;
                bf[0] += -x2 * (grad_trial[0] * (C2 * f[3] * x0 - grad_test[0] * (C1 + C2 * pow(f[3], 2))) +
                                grad_trial[1] * (C2 * f[2] * x1 - grad_test[1] * (C1 + C2 * pow(f[2], 2))));
                bf[1] += -dx * (grad_trial[0] * (grad_test[1] * x7 + x4) + grad_trial[1] * (-grad_test[0] * x9 + x8));
                bf[3] += -dx * (grad_trial[0] * (-grad_test[1] * x9 + x4) + grad_trial[1] * (grad_test[0] * x7 + x8));
                bf[4] += -x2 * (grad_trial[0] * (C2 * f[1] * x10 - grad_test[0] * (C1 + C2 * pow(f[1], 2))) +
                                grad_trial[1] * (C2 * f[0] * x11 - grad_test[1] * (C1 + C2 * pow(f[0], 2))));
                bf[8] += 0;
                bf[6] += x12 * (-f[2] * grad_trial[1] + f[3] * grad_trial[0]);
                bf[7] += x12 * (f[0] * grad_trial[1] - f[1] * grad_trial[0]);
                bf[2] += x13 * (-x0 + x1);
                bf[5] += x13 * (x10 - x11);
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T p,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T fun_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = grad_test[0] * p;
                T x1 = 2 * C1;
                T x2 = grad_test[0] * x1;
                T x3 = grad_test[1] * x1;
                T x4 = grad_test[1] * p;
                T x5 = f[0] * f[3];
                T x6 = 2 * C2;
                T x7 = grad_test[1] * x6;
                T x8 = f[1] * f[2];
                T x9 = grad_test[0] * x6;
                lf[0] += dx * (f[0] * pow(f[3], 2) * x9 + f[0] * x2 + f[1] * pow(f[2], 2) * x7 + f[1] * x3 - f[2] * x4 -
                               f[2] * x5 * x7 + f[3] * x0 - f[3] * x8 * x9);
                lf[1] += dx * (pow(f[0], 2) * f[3] * x7 + f[0] * x4 - f[0] * x7 * x8 + pow(f[1], 2) * f[2] * x9 -
                               f[1] * x0 - f[1] * x5 * x9 + f[2] * x2 + f[3] * x3);
                lf[2] += dx * fun_test * (x5 - x8 - 1);
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T p, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = pow(f[0], 2) + pow(f[2], 2);
                T x1 = pow(f[1], 2) + pow(f[3], 2);
                T x2 = x0 + x1;
                e += dx * (C1 * (x2 - 2) +
                           C2 * (-1.0 / 2.0 * pow(x0, 2) - 1.0 / 2.0 * pow(x1, 2) + (1.0 / 2.0) * pow(x2, 2) -
                                 pow(f[0] * f[1] + f[2] * f[3], 2) - 2) +
                           p * (f[0] * f[3] - f[1] * f[2] - 1));
            }

            UTOPIA_FUNCTION void eval(const T *UTOPIA_RESTRICT f,
                                      const T p,
                                      const T *grad_test,
                                      const T *grad_trial,
                                      const T fun_test,
                                      const T fun_trial,
                                      const T dx,
                                      T &e,
                                      T *UTOPIA_RESTRICT lf,
                                      T *UTOPIA_RESTRICT bf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[0] * f[3];
                T x1 = f[1] * f[2];
                T x2 = -x1;
                T x3 = x0 + x2 - 1;
                T x4 = pow(f[0], 2);
                T x5 = pow(f[2], 2);
                T x6 = x4 + x5;
                T x7 = pow(f[1], 2);
                T x8 = pow(f[3], 2);
                T x9 = x7 + x8;
                T x10 = x6 + x9;
                T x11 = f[0] * f[1];
                T x12 = f[2] * f[3];
                T x13 = f[3] * grad_test[0];
                T x14 = 2 * C1;
                T x15 = grad_test[0] * x14;
                T x16 = grad_test[1] * x14;
                T x17 = f[2] * grad_test[1];
                T x18 = f[0] * grad_test[1];
                T x19 = 2 * C2;
                T x20 = x12 * x19;
                T x21 = f[1] * grad_test[0];
                T x22 = C2 * x8;
                T x23 = 2 * grad_test[0];
                T x24 = C2 * x5;
                T x25 = 2 * grad_test[1];
                T x26 = x17 * x19;
                T x27 = x13 * x19;
                T x28 = C2 * x4;
                T x29 = C2 * x7;
                T x30 = dx * fun_test;
                T x31 = C2 * x12;
                T x32 = 2 * dx;
                T x33 = f[1] * x27;
                T x34 = p + x19 * (x0 - 2 * x1);
                T x35 = f[0] * x26;
                T x36 = p + x19 * (2 * x0 + x2);
                T x37 = C2 * x11;
                T x38 = dx * fun_trial;
                e += dx * (C1 * (x10 - 2) +
                           C2 * ((1.0 / 2.0) * pow(x10, 2) - 1.0 / 2.0 * pow(x6, 2) - 1.0 / 2.0 * pow(x9, 2) -
                                 pow(x11 + x12, 2) - 2) +
                           p * x3);
                lf[0] += dx * (f[0] * x15 + f[0] * x22 * x23 + f[1] * x16 + f[1] * x24 * x25 + p * x13 - p * x17 -
                               x18 * x20 - x20 * x21);
                lf[1] += dx * (f[2] * x15 + f[2] * x23 * x29 + f[3] * x16 + f[3] * x25 * x28 + p * x18 - p * x21 -
                               x11 * x26 - x11 * x27);
                lf[2] += x3 * x30;
                bf[0] += -x32 * (grad_trial[0] * (-grad_test[0] * (C1 + x22) + grad_test[1] * x31) +
                                 grad_trial[1] * (grad_test[0] * x31 - grad_test[1] * (C1 + x24)));
                bf[1] +=
                    -dx * (grad_trial[0] * (grad_test[1] * x34 + x33) + grad_trial[1] * (-grad_test[0] * x36 + x35));
                bf[3] +=
                    -dx * (grad_trial[0] * (-grad_test[1] * x36 + x33) + grad_trial[1] * (grad_test[0] * x34 + x35));
                bf[4] += -x32 * (grad_trial[0] * (-grad_test[0] * (C1 + x29) + grad_test[1] * x37) +
                                 grad_trial[1] * (grad_test[0] * x37 - grad_test[1] * (C1 + x28)));
                bf[8] += 0;
                bf[6] += x30 * (-f[2] * grad_trial[1] + f[3] * grad_trial[0]);
                bf[7] += x30 * (f[0] * grad_trial[1] - f[1] * grad_trial[0]);
                bf[2] += x38 * (x13 - x17);
                bf[5] += x38 * (x18 - x21);
            }

            T C1{1.0};
            T C2{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_2_IMPL_hpp
