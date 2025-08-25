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

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "IncompressibleMooneyRivlin_2"; }

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    in.get("C1", C1);
                    in.get("C2", C2);
                    in.get("K", K);
                }

                T C1{0.083};
                T C2{0.083};
                T K{166.67};
            };

            IncompressibleMooneyRivlin(const Params &params) {
                C1 = params.C1;
                C2 = params.C2;
                K = params.K;
            }

            UTOPIA_FUNCTION void hessian(const T *UTOPIA_RESTRICT f,
                                         const T *grad_test,
                                         const T *grad_trial,
                                         const T dx,
                                         T *UTOPIA_RESTRICT bf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = pow(f[3], 2);
                T x1 = pow(f[0], 2);
                T x2 = pow(f[1], 2);
                T x3 = x1 + x2;
                T x4 = pow(f[2], 2);
                T x5 = x0 + x4;
                T x6 = x3 + x5;
                T x7 = f[0] * f[3];
                T x8 = f[1] * f[2];
                T x9 = x7 - x8;
                T x10 = pow(x9, -8.0 / 3.0);
                T x11 = (10.0 / 9.0) * x10 * x6;
                T x12 = 2 / pow(x9, 2.0 / 3.0);
                T x13 = pow(x9, -5.0 / 3.0);
                T x14 = (8.0 / 3.0) * x13;
                T x15 = x12 - x14 * x7;
                T x16 = pow(x9, -4.0 / 3.0);
                T x17 = 2 * x16;
                T x18 = 2 * f[0];
                T x19 = f[0] * f[2];
                T x20 = f[1] * f[3];
                T x21 = x19 + x20;
                T x22 = 2 * x21;
                T x23 = 2 * f[0] * x6 - f[2] * x22 - x18 * x3;
                T x24 = pow(x9, -7.0 / 3.0);
                T x25 = (8.0 / 3.0) * x24;
                T x26 = -pow(x21, 2) - 1.0 / 2.0 * pow(x3, 2) - 1.0 / 2.0 * pow(x5, 2) + (1.0 / 2.0) * pow(x6, 2);
                T x27 = pow(x9, -10.0 / 3.0);
                T x28 = (28.0 / 9.0) * x26 * x27;
                T x29 = f[2] * f[3];
                T x30 = (4.0 / 3.0) * x13;
                T x31 = x19 * x30;
                T x32 = x20 * x30;
                T x33 = -2 * f[1] * x3 + 2 * f[1] * x6 - f[3] * x22;
                T x34 = (4.0 / 3.0) * x24;
                T x35 = f[3] * x34;
                T x36 = C1 * (-x11 * x29 + x31 - x32) +
                        C2 * ((4.0 / 3.0) * f[2] * x23 * x24 - x17 * x29 - x28 * x29 - x33 * x35) - K * x29;
                T x37 = x12 + x14 * x8;
                T x38 = f[0] * f[1];
                T x39 = x30 * x38;
                T x40 = x29 * x30;
                T x41 = -2 * f[2] * x5 + 2 * f[2] * x6 - x18 * x21;
                T x42 = grad_trial[0] *
                        (C1 * (-x11 * x20 + x39 - x40) +
                         C2 * ((4.0 / 3.0) * f[1] * x23 * x24 - x17 * x20 - x20 * x28 - x35 * x41) - K * x20);
                T x43 = K * (f[0] * f[3] - x8 - 1);
                T x44 = (2.0 / 3.0) * x13 * x6;
                T x45 = -f[1] * x22 - 2 * f[3] * x5 + 2 * f[3] * x6;
                T x46 = x26 * x34;
                T x47 = C1 * ((10.0 / 9.0) * f[0] * f[3] * x10 * x6 - x0 * x30 - x1 * x30 - x44) +
                        C2 * ((28.0 / 9.0) * f[0] * f[3] * x26 * x27 - f[0] * x23 * x34 + x16 * (4 * x7 - 2 * x8) -
                              x35 * x45 - x46) +
                        K * x7 + x43;
                T x48 = x33 * x34;
                T x49 = grad_trial[1] *
                        (C1 * (-x11 * x19 - x39 + x40) +
                         C2 * (-f[0] * x48 + (4.0 / 3.0) * f[2] * x24 * x45 - x17 * x19 - x19 * x28) - K * x19);
                T x50 = C1 * (x11 * x8 + x2 * x30 + x30 * x4 + x44) +
                        C2 * (f[1] * x48 + f[2] * x34 * x41 + x16 * (4 * f[1] * f[2] - 2 * x7) + x28 * x8 + x46) +
                        K * x8 - x43;
                T x51 = C1 * (-x11 * x38 - x31 + x32) +
                        C2 * (-f[0] * x34 * x41 + (4.0 / 3.0) * f[1] * x24 * x45 - x17 * x38 - x28 * x38) - K * x38;
                bf[0] +=
                    dx * (grad_test[0] * (grad_trial[0] * (C1 * (x0 * x11 + x15) +
                                                           C2 * (-f[3] * x23 * x25 + x0 * x17 + x0 * x28) + K * x0) +
                                          grad_trial[1] * x36) +
                          grad_test[1] * (grad_trial[0] * x36 +
                                          grad_trial[1] * (C1 * (x11 * x4 + x37) +
                                                           C2 * (f[2] * x25 * x33 + x17 * x4 + x28 * x4) + K * x4)));
                bf[1] += dx * (grad_test[0] * (grad_trial[1] * x47 + x42) + grad_test[1] * (grad_trial[0] * x50 + x49));
                bf[2] += dx * (grad_test[0] * (grad_trial[1] * x50 + x42) + grad_test[1] * (grad_trial[0] * x47 + x49));
                bf[3] +=
                    dx * (grad_test[0] * (grad_trial[0] * (C1 * (x11 * x2 + x37) +
                                                           C2 * (f[1] * x25 * x41 + x17 * x2 + x2 * x28) + K * x2) +
                                          grad_trial[1] * x51) +
                          grad_test[1] * (grad_trial[0] * x51 +
                                          grad_trial[1] * (C1 * (x1 * x11 + x15) +
                                                           C2 * (-f[0] * x25 * x45 + x1 * x17 + x1 * x28) + K * x1)));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[1] * f[2];
                T x1 = K * (f[0] * f[3] - x0 - 1);
                T x2 = f[0] * f[3] - x0;
                T x3 = pow(x2, -2.0 / 3.0);
                T x4 = 2 * f[0];
                T x5 = pow(f[0], 2) + pow(f[1], 2);
                T x6 = pow(f[2], 2) + pow(f[3], 2);
                T x7 = x5 + x6;
                T x8 = (2.0 / 3.0) * x7 / pow(x2, 5.0 / 3.0);
                T x9 = pow(x2, -4.0 / 3.0);
                T x10 = f[0] * f[2] + f[1] * f[3];
                T x11 = 2 * x10;
                T x12 = (4.0 / 3.0) *
                        (-pow(x10, 2) - 1.0 / 2.0 * pow(x5, 2) - 1.0 / 2.0 * pow(x6, 2) + (1.0 / 2.0) * pow(x7, 2)) /
                        pow(x2, 7.0 / 3.0);
                T x13 = 2 * f[1];
                T x14 = 2 * f[2];
                lf[0] +=
                    dx *
                    (grad_test[0] * (C1 * (-f[3] * x8 + x3 * x4) +
                                     C2 * (-f[3] * x12 + x9 * (2 * f[0] * x7 - f[2] * x11 - x4 * x5)) + f[3] * x1) +
                     grad_test[1] * (C1 * (f[2] * x8 + x13 * x3) +
                                     C2 * (f[2] * x12 + x9 * (2 * f[1] * x7 - f[3] * x11 - x13 * x5)) - f[2] * x1));
                lf[1] +=
                    dx * (grad_test[0] * (C1 * (f[1] * x8 + x14 * x3) +
                                          C2 * (f[1] * x12 + x9 * (2 * f[2] * x7 - x10 * x4 - x14 * x6)) - f[1] * x1) +
                          grad_test[1] *
                              (C1 * (-f[0] * x8 + 2 * f[3] * x3) +
                               C2 * (-f[0] * x12 + x9 * (-f[1] * x11 - 2 * f[3] * x6 + 2 * f[3] * x7)) + f[0] * x1));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[1] * f[2];
                T x1 = f[0] * f[3] - x0;
                T x2 = pow(f[0], 2) + pow(f[1], 2);
                T x3 = pow(f[2], 2) + pow(f[3], 2);
                T x4 = x2 + x3;
                e += dx * (C1 * (-2 + x4 / pow(x1, 2.0 / 3.0)) +
                           C2 * (-1 + (-1.0 / 2.0 * pow(x2, 2) - 1.0 / 2.0 * pow(x3, 2) + (1.0 / 2.0) * pow(x4, 2) -
                                       pow(f[0] * f[2] + f[1] * f[3], 2)) /
                                          pow(x1, 4.0 / 3.0)) +
                           (1.0 / 2.0) * K * pow(f[0] * f[3] - x0 - 1, 2));
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
                T x0 = f[1] * f[2];
                T x1 = f[0] * f[3] - x0 - 1;
                T x2 = f[0] * f[3];
                T x3 = -x0 + x2;
                T x4 = pow(x3, -2.0 / 3.0);
                T x5 = pow(f[0], 2);
                T x6 = pow(f[1], 2);
                T x7 = x5 + x6;
                T x8 = pow(f[2], 2);
                T x9 = pow(f[3], 2);
                T x10 = x8 + x9;
                T x11 = x10 + x7;
                T x12 = pow(x3, -4.0 / 3.0);
                T x13 = f[0] * f[2];
                T x14 = f[1] * f[3];
                T x15 = x13 + x14;
                T x16 = -1.0 / 2.0 * pow(x10, 2) + (1.0 / 2.0) * pow(x11, 2) - pow(x15, 2) - 1.0 / 2.0 * pow(x7, 2);
                T x17 = K * x1;
                T x18 = 2 * x4;
                T x19 = pow(x3, -5.0 / 3.0);
                T x20 = (2.0 / 3.0) * x11 * x19;
                T x21 = 2 * f[0];
                T x22 = 2 * x15;
                T x23 = 2 * f[0] * x11 - f[2] * x22 - x21 * x7;
                T x24 = pow(x3, -7.0 / 3.0);
                T x25 = (4.0 / 3.0) * x24;
                T x26 = x16 * x25;
                T x27 = 2 * f[1] * x11 - 2 * f[1] * x7 - f[3] * x22;
                T x28 = -f[1] * x22 - 2 * f[3] * x10 + 2 * f[3] * x11;
                T x29 = 2 * f[2];
                T x30 = 2 * f[2] * x11 - x10 * x29 - x15 * x21;
                T x31 = pow(x3, -8.0 / 3.0);
                T x32 = (10.0 / 9.0) * x11 * x31;
                T x33 = (8.0 / 3.0) * x19;
                T x34 = x18 - x2 * x33;
                T x35 = 2 * x12;
                T x36 = (8.0 / 3.0) * x24;
                T x37 = pow(x3, -10.0 / 3.0);
                T x38 = (28.0 / 9.0) * x16 * x37;
                T x39 = f[2] * f[3];
                T x40 = (4.0 / 3.0) * x19;
                T x41 = x13 * x40;
                T x42 = x14 * x40;
                T x43 = f[3] * x25;
                T x44 = C1 * (-x32 * x39 + x41 - x42) +
                        C2 * ((4.0 / 3.0) * f[2] * x23 * x24 - f[3] * x12 * x29 - x27 * x43 - x38 * x39) - K * x39;
                T x45 = x0 * x33 + x18;
                T x46 = f[0] * f[1];
                T x47 = x40 * x46;
                T x48 = x39 * x40;
                T x49 = grad_trial[0] *
                        (C1 * (-x14 * x32 + x47 - x48) +
                         C2 * ((4.0 / 3.0) * f[1] * x23 * x24 - x14 * x35 - x14 * x38 - x30 * x43) - K * x14);
                T x50 = C1 * ((10.0 / 9.0) * f[0] * f[3] * x11 * x31 - x20 - x40 * x5 - x40 * x9) +
                        C2 * ((28.0 / 9.0) * f[0] * f[3] * x16 * x37 - f[0] * x23 * x25 + x12 * (-2 * x0 + 4 * x2) -
                              x26 - x28 * x43) +
                        K * x2 + x17;
                T x51 = x25 * x27;
                T x52 = grad_trial[1] *
                        (C1 * (-x13 * x32 - x47 + x48) +
                         C2 * (-f[0] * x51 + (4.0 / 3.0) * f[2] * x24 * x28 - x13 * x35 - x13 * x38) - K * x13);
                T x53 = C1 * (x0 * x32 + x20 + x40 * x6 + x40 * x8) +
                        C2 * (f[1] * x51 + f[2] * x25 * x30 + x0 * x38 + x12 * (4 * f[1] * f[2] - 2 * x2) + x26) +
                        K * x0 - x17;
                T x54 = C1 * (-x32 * x46 - x41 + x42) +
                        C2 * (-f[0] * x25 * x30 - f[1] * x12 * x21 + (4.0 / 3.0) * f[1] * x24 * x28 - x38 * x46) -
                        K * x46;
                e += dx * (C1 * (x11 * x4 - 2) + C2 * (x12 * x16 - 1) + (1.0 / 2.0) * K * pow(x1, 2));
                lf[0] +=
                    dx *
                    (grad_test[0] * (C1 * (f[0] * x18 - f[3] * x20) + C2 * (-f[3] * x26 + x12 * x23) + f[3] * x17) +
                     grad_test[1] * (C1 * (f[1] * x18 + f[2] * x20) + C2 * (f[2] * x26 + x12 * x27) - f[2] * x17));
                lf[1] +=
                    dx *
                    (grad_test[0] * (C1 * (f[1] * x20 + f[2] * x18) + C2 * (f[1] * x26 + x12 * x30) - f[1] * x17) +
                     grad_test[1] * (C1 * (-f[0] * x20 + 2 * f[3] * x4) + C2 * (-f[0] * x26 + x12 * x28) + f[0] * x17));
                bf[0] +=
                    dx * (grad_test[0] * (grad_trial[0] * (C1 * (x32 * x9 + x34) +
                                                           C2 * (-f[3] * x23 * x36 + x35 * x9 + x38 * x9) + K * x9) +
                                          grad_trial[1] * x44) +
                          grad_test[1] * (grad_trial[0] * x44 +
                                          grad_trial[1] * (C1 * (x32 * x8 + x45) +
                                                           C2 * (f[2] * x27 * x36 + x35 * x8 + x38 * x8) + K * x8)));
                bf[1] += dx * (grad_test[0] * (grad_trial[1] * x50 + x49) + grad_test[1] * (grad_trial[0] * x53 + x52));
                bf[2] += dx * (grad_test[0] * (grad_trial[1] * x53 + x49) + grad_test[1] * (grad_trial[0] * x50 + x52));
                bf[3] +=
                    dx * (grad_test[0] * (grad_trial[0] * (C1 * (x32 * x6 + x45) +
                                                           C2 * (f[1] * x30 * x36 + x35 * x6 + x38 * x6) + K * x6) +
                                          grad_trial[1] * x54) +
                          grad_test[1] * (grad_trial[0] * x54 +
                                          grad_trial[1] * (C1 * (x32 * x5 + x34) +
                                                           C2 * (-f[0] * x28 * x36 + x35 * x5 + x38 * x5) + K * x5)));
            }

            UTOPIA_FUNCTION void apply(const T *UTOPIA_RESTRICT f,
                                       const T *grad_test,
                                       const T *disp_grad,
                                       const T dx,
                                       T *UTOPIA_RESTRICT res) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = pow(f[3], 2);
                T x1 = pow(f[0], 2);
                T x2 = pow(f[1], 2);
                T x3 = x1 + x2;
                T x4 = pow(f[2], 2);
                T x5 = x0 + x4;
                T x6 = x3 + x5;
                T x7 = f[0] * f[3];
                T x8 = f[1] * f[2];
                T x9 = x7 - x8;
                T x10 = pow(x9, -8.0 / 3.0);
                T x11 = (10.0 / 9.0) * x10 * x6;
                T x12 = 2 / pow(x9, 2.0 / 3.0);
                T x13 = pow(x9, -5.0 / 3.0);
                T x14 = (8.0 / 3.0) * x13;
                T x15 = x12 - x14 * x7;
                T x16 = pow(x9, -4.0 / 3.0);
                T x17 = 2 * x16;
                T x18 = 2 * f[0];
                T x19 = f[0] * f[2];
                T x20 = f[1] * f[3];
                T x21 = x19 + x20;
                T x22 = 2 * x21;
                T x23 = 2 * f[0] * x6 - f[2] * x22 - x18 * x3;
                T x24 = pow(x9, -7.0 / 3.0);
                T x25 = (8.0 / 3.0) * x24;
                T x26 = -pow(x21, 2) - 1.0 / 2.0 * pow(x3, 2) - 1.0 / 2.0 * pow(x5, 2) + (1.0 / 2.0) * pow(x6, 2);
                T x27 = pow(x9, -10.0 / 3.0);
                T x28 = (28.0 / 9.0) * x26 * x27;
                T x29 = f[2] * f[3];
                T x30 = (4.0 / 3.0) * x13;
                T x31 = x19 * x30;
                T x32 = x20 * x30;
                T x33 = -2 * f[1] * x3 + 2 * f[1] * x6 - f[3] * x22;
                T x34 = (4.0 / 3.0) * x24;
                T x35 = f[3] * x34;
                T x36 = C1 * (-x11 * x29 + x31 - x32) +
                        C2 * ((4.0 / 3.0) * f[2] * x23 * x24 - x17 * x29 - x28 * x29 - x33 * x35) - K * x29;
                T x37 = f[0] * f[1];
                T x38 = x30 * x37;
                T x39 = x29 * x30;
                T x40 = -2 * f[2] * x5 + 2 * f[2] * x6 - x18 * x21;
                T x41 = C1 * (-x11 * x20 + x38 - x39) +
                        C2 * ((4.0 / 3.0) * f[1] * x23 * x24 - x17 * x20 - x20 * x28 - x35 * x40) - K * x20;
                T x42 = K * (f[0] * f[3] - x8 - 1);
                T x43 = (2.0 / 3.0) * x13 * x6;
                T x44 = -f[1] * x22 - 2 * f[3] * x5 + 2 * f[3] * x6;
                T x45 = x26 * x34;
                T x46 = C1 * ((10.0 / 9.0) * f[0] * f[3] * x10 * x6 - x0 * x30 - x1 * x30 - x43) +
                        C2 * ((28.0 / 9.0) * f[0] * f[3] * x26 * x27 - f[0] * x23 * x34 + x16 * (4 * x7 - 2 * x8) -
                              x35 * x44 - x45) +
                        K * x7 + x42;
                T x47 = x12 + x14 * x8;
                T x48 = x33 * x34;
                T x49 = C1 * (-x11 * x19 - x38 + x39) +
                        C2 * (-f[0] * x48 + (4.0 / 3.0) * f[2] * x24 * x44 - x17 * x19 - x19 * x28) - K * x19;
                T x50 = C1 * (x11 * x8 + x2 * x30 + x30 * x4 + x43) +
                        C2 * (f[1] * x48 + f[2] * x34 * x40 + x16 * (4 * f[1] * f[2] - 2 * x7) + x28 * x8 + x45) +
                        K * x8 - x42;
                T x51 = C1 * (-x11 * x37 - x31 + x32) +
                        C2 * (-f[0] * x34 * x40 + (4.0 / 3.0) * f[1] * x24 * x44 - x17 * x37 - x28 * x37) - K * x37;
                res[0] +=
                    dx * (grad_test[0] * (disp_grad[0] * (C1 * (x0 * x11 + x15) +
                                                          C2 * (-f[3] * x23 * x25 + x0 * x17 + x0 * x28) + K * x0) +
                                          disp_grad[1] * x36 + disp_grad[2] * x41 + disp_grad[3] * x46) +
                          grad_test[1] * (disp_grad[0] * x36 +
                                          disp_grad[1] * (C1 * (x11 * x4 + x47) +
                                                          C2 * (f[2] * x25 * x33 + x17 * x4 + x28 * x4) + K * x4) +
                                          disp_grad[2] * x50 + disp_grad[3] * x49));
                res[1] +=
                    dx * (grad_test[0] * (disp_grad[0] * x41 + disp_grad[1] * x50 +
                                          disp_grad[2] * (C1 * (x11 * x2 + x47) +
                                                          C2 * (f[1] * x25 * x40 + x17 * x2 + x2 * x28) + K * x2) +
                                          disp_grad[3] * x51) +
                          grad_test[1] * (disp_grad[0] * x46 + disp_grad[1] * x49 + disp_grad[2] * x51 +
                                          disp_grad[3] * (C1 * (x1 * x11 + x15) +
                                                          C2 * (-f[0] * x25 * x44 + x1 * x17 + x1 * x28) + K * x1)));
            }

            T C1{0.083};
            T C2{0.083};
            T K{166.67};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_IncompressibleMooneyRivlin_2_IMPL_hpp
