#ifndef UTOPIA_TPL_HYPERELASTICITY_MooneyRivlin_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_MooneyRivlin_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_MooneyRivlin.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of MooneyRivlin for dimension 2
         */
        template <typename T>
        class MooneyRivlin<T, 2> {
        public:
            static constexpr int Dim = 2;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "MooneyRivlin_2"; }

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

            MooneyRivlin(const Params &params) {
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
                T x1 = f[0] * f[3];
                T x2 = f[1] * f[2];
                T x3 = x1 - x2;
                T x4 = K / pow(x3, 2);
                T x5 = x0 * x4;
                T x6 = log(x3);
                T x7 = pow(f[0], 2);
                T x8 = pow(f[1], 2);
                T x9 = x7 + x8;
                T x10 = pow(f[2], 2);
                T x11 = x0 + x10;
                T x12 = x11 + x9;
                T x13 = pow(x3, -8.0 / 3.0);
                T x14 = (10.0 / 9.0) * x12 * x13;
                T x15 = 2 / pow(x3, 2.0 / 3.0);
                T x16 = pow(x3, -5.0 / 3.0);
                T x17 = (8.0 / 3.0) * x16;
                T x18 = -x1 * x17 + x15;
                T x19 = pow(x3, -4.0 / 3.0);
                T x20 = 2 * x19;
                T x21 = 2 * f[0];
                T x22 = f[0] * f[2];
                T x23 = f[1] * f[3];
                T x24 = x22 + x23;
                T x25 = 2 * x24;
                T x26 = 2 * f[0] * x12 - f[2] * x25 - x21 * x9;
                T x27 = pow(x3, -7.0 / 3.0);
                T x28 = (8.0 / 3.0) * x27;
                T x29 = -1.0 / 2.0 * pow(x11, 2) + (1.0 / 2.0) * pow(x12, 2) - pow(x24, 2) - 1.0 / 2.0 * pow(x9, 2);
                T x30 = pow(x3, -10.0 / 3.0);
                T x31 = (28.0 / 9.0) * x29 * x30;
                T x32 = f[2] * f[3];
                T x33 = x32 * x4;
                T x34 = (4.0 / 3.0) * x16;
                T x35 = x22 * x34;
                T x36 = x23 * x34;
                T x37 = 2 * f[1] * x12 - 2 * f[1] * x9 - f[3] * x25;
                T x38 = (4.0 / 3.0) * x27;
                T x39 = f[3] * x38;
                T x40 = C1 * (-x14 * x32 + x35 - x36) +
                        C2 * ((4.0 / 3.0) * f[2] * x26 * x27 - x20 * x32 - x31 * x32 - x37 * x39) + x33 * x6 - x33;
                T x41 = x10 * x4;
                T x42 = x15 + x17 * x2;
                T x43 = x23 * x4;
                T x44 = f[0] * f[1];
                T x45 = x34 * x44;
                T x46 = x32 * x34;
                T x47 = -2 * f[2] * x11 + 2 * f[2] * x12 - x21 * x24;
                T x48 = grad_trial[0] *
                        (C1 * (-x14 * x23 + x45 - x46) +
                         C2 * ((4.0 / 3.0) * f[1] * x26 * x27 - x20 * x23 - x23 * x31 - x39 * x47) + x43 * x6 - x43);
                T x49 = x1 * x4;
                T x50 = K * x6 / x3;
                T x51 = (2.0 / 3.0) * x12 * x16;
                T x52 = -f[1] * x25 - 2 * f[3] * x11 + 2 * f[3] * x12;
                T x53 = x29 * x38;
                T x54 = C1 * ((10.0 / 9.0) * f[0] * f[3] * x12 * x13 - x0 * x34 - x34 * x7 - x51) +
                        C2 * ((28.0 / 9.0) * f[0] * f[3] * x29 * x30 - f[0] * x26 * x38 + x19 * (4 * x1 - 2 * x2) -
                              x39 * x52 - x53) -
                        x49 * x6 + x49 + x50;
                T x55 = x22 * x4;
                T x56 = x37 * x38;
                T x57 = grad_trial[1] *
                        (C1 * (-x14 * x22 - x45 + x46) +
                         C2 * (-f[0] * x56 + (4.0 / 3.0) * f[2] * x27 * x52 - x20 * x22 - x22 * x31) + x55 * x6 - x55);
                T x58 = x2 * x4;
                T x59 = C1 * (x10 * x34 + x14 * x2 + x34 * x8 + x51) +
                        C2 * (f[1] * x56 + f[2] * x38 * x47 + x19 * (4 * f[1] * f[2] - 2 * x1) + x2 * x31 + x53) - x50 -
                        x58 * x6 + x58;
                T x60 = x4 * x8;
                T x61 = x4 * x44;
                T x62 = C1 * (-x14 * x44 - x35 + x36) +
                        C2 * (-f[0] * x38 * x47 + (4.0 / 3.0) * f[1] * x27 * x52 - x20 * x44 - x31 * x44) + x6 * x61 -
                        x61;
                T x63 = x4 * x7;
                bf[0] +=
                    dx *
                    (grad_test[0] * (grad_trial[0] * (C1 * (x0 * x14 + x18) +
                                                      C2 * (-f[3] * x26 * x28 + x0 * x20 + x0 * x31) - x5 * x6 + x5) +
                                     grad_trial[1] * x40) +
                     grad_test[1] *
                         (grad_trial[0] * x40 +
                          grad_trial[1] * (C1 * (x10 * x14 + x42) + C2 * (f[2] * x28 * x37 + x10 * x20 + x10 * x31) -
                                           x41 * x6 + x41)));
                bf[1] += dx * (grad_test[0] * (grad_trial[1] * x54 + x48) + grad_test[1] * (grad_trial[0] * x59 + x57));
                bf[2] += dx * (grad_test[0] * (grad_trial[1] * x59 + x48) + grad_test[1] * (grad_trial[0] * x54 + x57));
                bf[3] += dx * (grad_test[0] *
                                   (grad_trial[0] * (C1 * (x14 * x8 + x42) +
                                                     C2 * (f[1] * x28 * x47 + x20 * x8 + x31 * x8) - x6 * x60 + x60) +
                                    grad_trial[1] * x62) +
                               grad_test[1] *
                                   (grad_trial[0] * x62 +
                                    grad_trial[1] * (C1 * (x14 * x7 + x18) +
                                                     C2 * (-f[0] * x28 * x52 + x20 * x7 + x31 * x7) - x6 * x63 + x63)));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[0] * f[3] - f[1] * f[2];
                T x1 = K * log(x0) / x0;
                T x2 = pow(x0, -2.0 / 3.0);
                T x3 = 2 * f[0];
                T x4 = pow(f[0], 2) + pow(f[1], 2);
                T x5 = pow(f[2], 2) + pow(f[3], 2);
                T x6 = x4 + x5;
                T x7 = (2.0 / 3.0) * x6 / pow(x0, 5.0 / 3.0);
                T x8 = pow(x0, -4.0 / 3.0);
                T x9 = f[0] * f[2] + f[1] * f[3];
                T x10 = 2 * x9;
                T x11 = (4.0 / 3.0) *
                        (-1.0 / 2.0 * pow(x4, 2) - 1.0 / 2.0 * pow(x5, 2) + (1.0 / 2.0) * pow(x6, 2) - pow(x9, 2)) /
                        pow(x0, 7.0 / 3.0);
                T x12 = 2 * f[1];
                T x13 = 2 * f[2];
                lf[0] +=
                    dx *
                    (grad_test[0] * (C1 * (-f[3] * x7 + x2 * x3) +
                                     C2 * (-f[3] * x11 + x8 * (2 * f[0] * x6 - f[2] * x10 - x3 * x4)) + f[3] * x1) +
                     grad_test[1] * (C1 * (f[2] * x7 + x12 * x2) +
                                     C2 * (f[2] * x11 + x8 * (2 * f[1] * x6 - f[3] * x10 - x12 * x4)) - f[2] * x1));
                lf[1] +=
                    dx * (grad_test[0] * (C1 * (f[1] * x7 + x13 * x2) +
                                          C2 * (f[1] * x11 + x8 * (2 * f[2] * x6 - x13 * x5 - x3 * x9)) - f[1] * x1) +
                          grad_test[1] *
                              (C1 * (-f[0] * x7 + 2 * f[3] * x2) +
                               C2 * (-f[0] * x11 + x8 * (-f[1] * x10 - 2 * f[3] * x5 + 2 * f[3] * x6)) + f[0] * x1));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[0] * f[3] - f[1] * f[2];
                T x1 = pow(f[0], 2) + pow(f[1], 2);
                T x2 = pow(f[2], 2) + pow(f[3], 2);
                T x3 = x1 + x2;
                e += dx * (C1 * (-2 + x3 / pow(x0, 2.0 / 3.0)) +
                           C2 * (-1 + (-1.0 / 2.0 * pow(x1, 2) - 1.0 / 2.0 * pow(x2, 2) + (1.0 / 2.0) * pow(x3, 2) -
                                       pow(f[0] * f[2] + f[1] * f[3], 2)) /
                                          pow(x0, 4.0 / 3.0)) +
                           (1.0 / 2.0) * K * pow(log(x0), 2));
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
                T x4 = pow(x2, -2.0 / 3.0);
                T x5 = pow(f[0], 2);
                T x6 = pow(f[1], 2);
                T x7 = x5 + x6;
                T x8 = pow(f[2], 2);
                T x9 = pow(f[3], 2);
                T x10 = x8 + x9;
                T x11 = x10 + x7;
                T x12 = pow(x2, -4.0 / 3.0);
                T x13 = f[0] * f[2];
                T x14 = f[1] * f[3];
                T x15 = x13 + x14;
                T x16 = -1.0 / 2.0 * pow(x10, 2) + (1.0 / 2.0) * pow(x11, 2) - pow(x15, 2) - 1.0 / 2.0 * pow(x7, 2);
                T x17 = K * x3 / x2;
                T x18 = 2 * x4;
                T x19 = pow(x2, -5.0 / 3.0);
                T x20 = (2.0 / 3.0) * x11 * x19;
                T x21 = 2 * f[0];
                T x22 = 2 * x15;
                T x23 = 2 * f[0] * x11 - f[2] * x22 - x21 * x7;
                T x24 = pow(x2, -7.0 / 3.0);
                T x25 = (4.0 / 3.0) * x24;
                T x26 = x16 * x25;
                T x27 = 2 * f[1] * x11 - 2 * f[1] * x7 - f[3] * x22;
                T x28 = -f[1] * x22 - 2 * f[3] * x10 + 2 * f[3] * x11;
                T x29 = 2 * f[2];
                T x30 = 2 * f[2] * x11 - x10 * x29 - x15 * x21;
                T x31 = K / pow(x2, 2);
                T x32 = x31 * x9;
                T x33 = pow(x2, -8.0 / 3.0);
                T x34 = (10.0 / 9.0) * x11 * x33;
                T x35 = (8.0 / 3.0) * x19;
                T x36 = -x0 * x35 + x18;
                T x37 = 2 * x12;
                T x38 = (8.0 / 3.0) * x24;
                T x39 = pow(x2, -10.0 / 3.0);
                T x40 = (28.0 / 9.0) * x16 * x39;
                T x41 = f[2] * f[3];
                T x42 = x31 * x41;
                T x43 = (4.0 / 3.0) * x19;
                T x44 = x13 * x43;
                T x45 = x14 * x43;
                T x46 = f[3] * x25;
                T x47 = C1 * (-x34 * x41 + x44 - x45) +
                        C2 * ((4.0 / 3.0) * f[2] * x23 * x24 - f[3] * x12 * x29 - x27 * x46 - x40 * x41) + x3 * x42 -
                        x42;
                T x48 = x31 * x8;
                T x49 = x1 * x35 + x18;
                T x50 = x14 * x31;
                T x51 = f[0] * f[1];
                T x52 = x43 * x51;
                T x53 = x41 * x43;
                T x54 = grad_trial[0] *
                        (C1 * (-x14 * x34 + x52 - x53) +
                         C2 * ((4.0 / 3.0) * f[1] * x23 * x24 - x14 * x37 - x14 * x40 - x30 * x46) + x3 * x50 - x50);
                T x55 = x0 * x31;
                T x56 = C1 * ((10.0 / 9.0) * f[0] * f[3] * x11 * x33 - x20 - x43 * x5 - x43 * x9) +
                        C2 * ((28.0 / 9.0) * f[0] * f[3] * x16 * x39 - f[0] * x23 * x25 + x12 * (4 * x0 - 2 * x1) -
                              x26 - x28 * x46) +
                        x17 - x3 * x55 + x55;
                T x57 = x13 * x31;
                T x58 = x25 * x27;
                T x59 = grad_trial[1] *
                        (C1 * (-x13 * x34 - x52 + x53) +
                         C2 * (-f[0] * x58 + (4.0 / 3.0) * f[2] * x24 * x28 - x13 * x37 - x13 * x40) + x3 * x57 - x57);
                T x60 = x1 * x31;
                T x61 = C1 * (x1 * x34 + x20 + x43 * x6 + x43 * x8) +
                        C2 * (f[1] * x58 + f[2] * x25 * x30 + x1 * x40 + x12 * (4 * f[1] * f[2] - 2 * x0) + x26) - x17 -
                        x3 * x60 + x60;
                T x62 = x31 * x6;
                T x63 = x31 * x51;
                T x64 = C1 * (-x34 * x51 - x44 + x45) +
                        C2 * (-f[0] * x25 * x30 - f[1] * x12 * x21 + (4.0 / 3.0) * f[1] * x24 * x28 - x40 * x51) +
                        x3 * x63 - x63;
                T x65 = x31 * x5;
                e += dx * (C1 * (x11 * x4 - 2) + C2 * (x12 * x16 - 1) + (1.0 / 2.0) * K * pow(x3, 2));
                lf[0] +=
                    dx *
                    (grad_test[0] * (C1 * (f[0] * x18 - f[3] * x20) + C2 * (-f[3] * x26 + x12 * x23) + f[3] * x17) +
                     grad_test[1] * (C1 * (f[1] * x18 + f[2] * x20) + C2 * (f[2] * x26 + x12 * x27) - f[2] * x17));
                lf[1] +=
                    dx *
                    (grad_test[0] * (C1 * (f[1] * x20 + f[2] * x18) + C2 * (f[1] * x26 + x12 * x30) - f[1] * x17) +
                     grad_test[1] * (C1 * (-f[0] * x20 + 2 * f[3] * x4) + C2 * (-f[0] * x26 + x12 * x28) + f[0] * x17));
                bf[0] +=
                    dx *
                    (grad_test[0] * (grad_trial[0] * (C1 * (x34 * x9 + x36) +
                                                      C2 * (-f[3] * x23 * x38 + x37 * x9 + x40 * x9) - x3 * x32 + x32) +
                                     grad_trial[1] * x47) +
                     grad_test[1] * (grad_trial[0] * x47 +
                                     grad_trial[1] * (C1 * (x34 * x8 + x49) +
                                                      C2 * (f[2] * x27 * x38 + x37 * x8 + x40 * x8) - x3 * x48 + x48)));
                bf[1] += dx * (grad_test[0] * (grad_trial[1] * x56 + x54) + grad_test[1] * (grad_trial[0] * x61 + x59));
                bf[2] += dx * (grad_test[0] * (grad_trial[1] * x61 + x54) + grad_test[1] * (grad_trial[0] * x56 + x59));
                bf[3] += dx * (grad_test[0] *
                                   (grad_trial[0] * (C1 * (x34 * x6 + x49) +
                                                     C2 * (f[1] * x30 * x38 + x37 * x6 + x40 * x6) - x3 * x62 + x62) +
                                    grad_trial[1] * x64) +
                               grad_test[1] *
                                   (grad_trial[0] * x64 +
                                    grad_trial[1] * (C1 * (x34 * x5 + x36) +
                                                     C2 * (-f[0] * x28 * x38 + x37 * x5 + x40 * x5) - x3 * x65 + x65)));
            }

            UTOPIA_FUNCTION void apply(const T *UTOPIA_RESTRICT f,
                                       const T *grad_test,
                                       const T *disp_grad,
                                       const T dx,
                                       T *UTOPIA_RESTRICT res) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = pow(f[3], 2);
                T x1 = f[0] * f[3];
                T x2 = f[1] * f[2];
                T x3 = x1 - x2;
                T x4 = K / pow(x3, 2);
                T x5 = x0 * x4;
                T x6 = log(x3);
                T x7 = pow(f[0], 2);
                T x8 = pow(f[1], 2);
                T x9 = x7 + x8;
                T x10 = pow(f[2], 2);
                T x11 = x0 + x10;
                T x12 = x11 + x9;
                T x13 = pow(x3, -8.0 / 3.0);
                T x14 = (10.0 / 9.0) * x12 * x13;
                T x15 = 2 / pow(x3, 2.0 / 3.0);
                T x16 = pow(x3, -5.0 / 3.0);
                T x17 = (8.0 / 3.0) * x16;
                T x18 = -x1 * x17 + x15;
                T x19 = pow(x3, -4.0 / 3.0);
                T x20 = 2 * x19;
                T x21 = 2 * f[0];
                T x22 = f[0] * f[2];
                T x23 = f[1] * f[3];
                T x24 = x22 + x23;
                T x25 = 2 * x24;
                T x26 = 2 * f[0] * x12 - f[2] * x25 - x21 * x9;
                T x27 = pow(x3, -7.0 / 3.0);
                T x28 = (8.0 / 3.0) * x27;
                T x29 = -1.0 / 2.0 * pow(x11, 2) + (1.0 / 2.0) * pow(x12, 2) - pow(x24, 2) - 1.0 / 2.0 * pow(x9, 2);
                T x30 = pow(x3, -10.0 / 3.0);
                T x31 = (28.0 / 9.0) * x29 * x30;
                T x32 = f[2] * f[3];
                T x33 = x32 * x4;
                T x34 = (4.0 / 3.0) * x16;
                T x35 = x22 * x34;
                T x36 = x23 * x34;
                T x37 = 2 * f[1] * x12 - 2 * f[1] * x9 - f[3] * x25;
                T x38 = (4.0 / 3.0) * x27;
                T x39 = f[3] * x38;
                T x40 = C1 * (-x14 * x32 + x35 - x36) +
                        C2 * ((4.0 / 3.0) * f[2] * x26 * x27 - x20 * x32 - x31 * x32 - x37 * x39) + x33 * x6 - x33;
                T x41 = x23 * x4;
                T x42 = f[0] * f[1];
                T x43 = x34 * x42;
                T x44 = x32 * x34;
                T x45 = -2 * f[2] * x11 + 2 * f[2] * x12 - x21 * x24;
                T x46 = C1 * (-x14 * x23 + x43 - x44) +
                        C2 * ((4.0 / 3.0) * f[1] * x26 * x27 - x20 * x23 - x23 * x31 - x39 * x45) + x41 * x6 - x41;
                T x47 = x1 * x4;
                T x48 = K * x6 / x3;
                T x49 = (2.0 / 3.0) * x12 * x16;
                T x50 = -f[1] * x25 - 2 * f[3] * x11 + 2 * f[3] * x12;
                T x51 = x29 * x38;
                T x52 = C1 * ((10.0 / 9.0) * f[0] * f[3] * x12 * x13 - x0 * x34 - x34 * x7 - x49) +
                        C2 * ((28.0 / 9.0) * f[0] * f[3] * x29 * x30 - f[0] * x26 * x38 + x19 * (4 * x1 - 2 * x2) -
                              x39 * x50 - x51) -
                        x47 * x6 + x47 + x48;
                T x53 = x10 * x4;
                T x54 = x15 + x17 * x2;
                T x55 = x22 * x4;
                T x56 = x37 * x38;
                T x57 = C1 * (-x14 * x22 - x43 + x44) +
                        C2 * (-f[0] * x56 + (4.0 / 3.0) * f[2] * x27 * x50 - x20 * x22 - x22 * x31) + x55 * x6 - x55;
                T x58 = x2 * x4;
                T x59 = C1 * (x10 * x34 + x14 * x2 + x34 * x8 + x49) +
                        C2 * (f[1] * x56 + f[2] * x38 * x45 + x19 * (4 * f[1] * f[2] - 2 * x1) + x2 * x31 + x51) - x48 -
                        x58 * x6 + x58;
                T x60 = x4 * x7;
                T x61 = x4 * x42;
                T x62 = C1 * (-x14 * x42 - x35 + x36) +
                        C2 * (-f[0] * x38 * x45 + (4.0 / 3.0) * f[1] * x27 * x50 - x20 * x42 - x31 * x42) + x6 * x61 -
                        x61;
                T x63 = x4 * x8;
                res[0] +=
                    dx *
                    (grad_test[0] * (disp_grad[0] * (C1 * (x0 * x14 + x18) +
                                                     C2 * (-f[3] * x26 * x28 + x0 * x20 + x0 * x31) - x5 * x6 + x5) +
                                     disp_grad[1] * x40 + disp_grad[2] * x46 + disp_grad[3] * x52) +
                     grad_test[1] * (disp_grad[0] * x40 +
                                     disp_grad[1] * (C1 * (x10 * x14 + x54) +
                                                     C2 * (f[2] * x28 * x37 + x10 * x20 + x10 * x31) - x53 * x6 + x53) +
                                     disp_grad[2] * x59 + disp_grad[3] * x57));
                res[1] +=
                    dx *
                    (grad_test[0] * (disp_grad[0] * x46 + disp_grad[1] * x59 +
                                     disp_grad[2] * (C1 * (x14 * x8 + x54) +
                                                     C2 * (f[1] * x28 * x45 + x20 * x8 + x31 * x8) - x6 * x63 + x63) +
                                     disp_grad[3] * x62) +
                     grad_test[1] * (disp_grad[0] * x52 + disp_grad[1] * x57 + disp_grad[2] * x62 +
                                     disp_grad[3] * (C1 * (x14 * x7 + x18) +
                                                     C2 * (-f[0] * x28 * x50 + x20 * x7 + x31 * x7) - x6 * x60 + x60)));
            }

            T C1{0.083};
            T C2{0.083};
            T K{166.67};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_MooneyRivlin_2_IMPL_hpp
