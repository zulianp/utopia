#ifndef UTOPIA_TPL_HYPERELASTICITY_Yeoh_2_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_Yeoh_2_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_Yeoh.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of Yeoh for dimension 2
         */
        template <typename T>
        class Yeoh<T, 2> {
        public:
            static constexpr int Dim = 2;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "Yeoh_2"; }

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    in.get("C0_0", C0_0);
                    in.get("C1_0", C1_0);
                    in.get("C0_1", C0_1);
                    in.get("C1_1", C1_1);
                }

                T C0_0{1.0};
                T C1_0{1.0};
                T C0_1{1.0};
                T C1_1{1.0};
            };

            Yeoh(const Params &params) {
                C0_0 = params.C0_0;
                C1_0 = params.C1_0;
                C0_1 = params.C0_1;
                C1_1 = params.C1_1;
            }

            UTOPIA_FUNCTION void hessian(const T *UTOPIA_RESTRICT f,
                                         const T *grad_test,
                                         const T *grad_trial,
                                         const T dx,
                                         T *UTOPIA_RESTRICT bf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = pow(f[3], 2);
                T x1 = 2 * C1_0;
                T x2 = f[0] * f[3];
                T x3 = f[1] * f[2];
                T x4 = x2 - x3;
                T x5 = x4 - 1;
                T x6 = 12 * C1_1 * pow(x5, 2);
                T x7 = pow(f[0], 2);
                T x8 = pow(f[1], 2);
                T x9 = pow(f[2], 2);
                T x10 = x0 + x7 + x8 + x9;
                T x11 = x10 / pow(x4, 8.0 / 3.0);
                T x12 = x0 * x11;
                T x13 = pow(x4, -2.0 / 3.0);
                T x14 = 2 * x13;
                T x15 = pow(x4, -5.0 / 3.0);
                T x16 = (8.0 / 3.0) * x15;
                T x17 = x14 - x16 * x2;
                T x18 = 4 * x13;
                T x19 = x10 * x15;
                T x20 = (4.0 / 3.0) * x19;
                T x21 = f[0] * x18 - f[3] * x20;
                T x22 = (2.0 / 3.0) * x19;
                T x23 = C0_1 * (f[0] * x14 - f[3] * x22);
                T x24 = (16.0 / 3.0) * x15;
                T x25 = x18 - x2 * x24;
                T x26 = C0_1 * (x10 * x13 - 2);
                T x27 = f[1] * x18 + f[2] * x20;
                T x28 = (4.0 / 3.0) * x15;
                T x29 = f[0] * f[2];
                T x30 = x28 * x29;
                T x31 = f[1] * f[3];
                T x32 = x28 * x31;
                T x33 = f[2] * f[3];
                T x34 = x11 * x33;
                T x35 = x16 * x29;
                T x36 = x16 * x31;
                T x37 =
                    C0_0 * (x30 - x32 - 10.0 / 9.0 * x34) - x1 * x33 + x26 * (-20.0 / 9.0 * x34 + x35 - x36) - x33 * x6;
                T x38 = C0_1 * (f[1] * x14 + f[2] * x22);
                T x39 = x11 * x9;
                T x40 = x14 + x16 * x3;
                T x41 = x18 + x24 * x3;
                T x42 = C0_1 * (f[1] * x22 + f[2] * x14);
                T x43 = f[0] * f[1];
                T x44 = x28 * x43;
                T x45 = x28 * x33;
                T x46 = x11 * x31;
                T x47 = x16 * x43;
                T x48 = x16 * x33;
                T x49 =
                    C0_0 * (x44 - x45 - 10.0 / 9.0 * x46) - x1 * x31 + x26 * (-20.0 / 9.0 * x46 + x47 - x48) - x31 * x6;
                T x50 = x11 * x3;
                T x51 = 4 * C1_1 * pow(x5, 3);
                T x52 = x1 * x5;
                T x53 = C0_0 * (x22 + x28 * x8 + x28 * x9 + (10.0 / 9.0) * x50) + x1 * x3 +
                        x26 * (x16 * x8 + x16 * x9 + x20 + (20.0 / 9.0) * x50) + x3 * x6 - x51 - x52;
                T x54 = C0_1 * (-f[0] * x22 + f[3] * x14);
                T x55 = x11 * x29;
                T x56 = C0_0 * (-x44 + x45 - 10.0 / 9.0 * x55) - x1 * x29 + x26 * (-x47 + x48 - 20.0 / 9.0 * x55) -
                        x29 * x6;
                T x57 = x11 * x2;
                T x58 = C0_0 * (-x0 * x28 - x22 - x28 * x7 + (10.0 / 9.0) * x57) + x1 * x2 + x2 * x6 +
                        x26 * (-x0 * x16 - x16 * x7 - x20 + (20.0 / 9.0) * x57) + x51 + x52;
                T x59 = f[1] * x20 + f[2] * x18;
                T x60 = -f[0] * x20 + f[3] * x18;
                T x61 = x11 * x8;
                T x62 = x11 * x43;
                T x63 = C0_0 * (-x30 + x32 - 10.0 / 9.0 * x62) - x1 * x43 + x26 * (-x35 + x36 - 20.0 / 9.0 * x62) -
                        x43 * x6;
                T x64 = x11 * x7;
                bf[0] +=
                    dx * (grad_trial[0] * (grad_test[0] * (C0_0 * ((10.0 / 9.0) * x12 + x17) + x0 * x1 + x0 * x6 +
                                                           x21 * x23 + x26 * ((20.0 / 9.0) * x12 + x25)) +
                                           grad_test[1] * (x23 * x27 + x37)) +
                          grad_trial[1] * (grad_test[0] * (x21 * x38 + x37) +
                                           grad_test[1] * (C0_0 * ((10.0 / 9.0) * x39 + x40) + x1 * x9 +
                                                           x26 * ((20.0 / 9.0) * x39 + x41) + x27 * x38 + x6 * x9)));
                bf[1] += dx * (grad_trial[0] * (grad_test[0] * (x21 * x42 + x49) + grad_test[1] * (x27 * x42 + x53)) +
                               grad_trial[1] * (grad_test[0] * (x21 * x54 + x58) + grad_test[1] * (x27 * x54 + x56)));
                bf[2] += dx * (grad_trial[0] * (grad_test[0] * (x23 * x59 + x49) + grad_test[1] * (x23 * x60 + x58)) +
                               grad_trial[1] * (grad_test[0] * (x38 * x59 + x53) + grad_test[1] * (x38 * x60 + x56)));
                bf[3] +=
                    dx * (grad_trial[0] * (grad_test[0] * (C0_0 * (x40 + (10.0 / 9.0) * x61) + x1 * x8 +
                                                           x26 * (x41 + (20.0 / 9.0) * x61) + x42 * x59 + x6 * x8) +
                                           grad_test[1] * (x42 * x60 + x63)) +
                          grad_trial[1] * (grad_test[0] * (x54 * x59 + x63) +
                                           grad_test[1] * (C0_0 * (x17 + (10.0 / 9.0) * x64) + x1 * x7 +
                                                           x26 * (x25 + (20.0 / 9.0) * x64) + x54 * x60 + x6 * x7)));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[0] * f[3] - f[1] * f[2];
                T x1 = x0 - 1;
                T x2 = 2 * C1_0 * x1;
                T x3 = 4 * C1_1 * pow(x1, 3);
                T x4 = pow(x0, -2.0 / 3.0);
                T x5 = f[0] * x4;
                T x6 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2);
                T x7 = x6 / pow(x0, 5.0 / 3.0);
                T x8 = f[3] * x7;
                T x9 = C0_1 * (x4 * x6 - 2);
                T x10 = f[1] * x4;
                T x11 = f[2] * x7;
                T x12 = f[2] * x4;
                T x13 = f[1] * x7;
                T x14 = f[3] * x4;
                T x15 = f[0] * x7;
                lf[0] += dx * (grad_test[0] * (C0_0 * (2 * x5 - 2.0 / 3.0 * x8) + f[3] * x2 + f[3] * x3 +
                                               x9 * (4 * x5 - 4.0 / 3.0 * x8)) +
                               grad_test[1] * (C0_0 * (2 * x10 + (2.0 / 3.0) * x11) - f[2] * x2 - f[2] * x3 +
                                               x9 * (4 * x10 + (4.0 / 3.0) * x11)));
                lf[1] += dx * (grad_test[0] * (C0_0 * (2 * x12 + (2.0 / 3.0) * x13) - f[1] * x2 - f[1] * x3 +
                                               x9 * (4 * x12 + (4.0 / 3.0) * x13)) +
                               grad_test[1] * (C0_0 * (2 * x14 - 2.0 / 3.0 * x15) + f[0] * x2 + f[0] * x3 +
                                               x9 * (4 * x14 - 4.0 / 3.0 * x15)));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[0] * f[3] - f[1] * f[2];
                T x1 = x0 - 1;
                T x2 = -2 + (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2)) / pow(x0, 2.0 / 3.0);
                e += dx * (C0_0 * x2 + C0_1 * pow(x2, 2) + C1_0 * pow(x1, 2) + C1_1 * pow(x1, 4));
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
                T x3 = x2 - 1;
                T x4 = pow(x3, 2);
                T x5 = pow(x2, -2.0 / 3.0);
                T x6 = pow(f[0], 2);
                T x7 = pow(f[1], 2);
                T x8 = pow(f[2], 2);
                T x9 = pow(f[3], 2);
                T x10 = x6 + x7 + x8 + x9;
                T x11 = x10 * x5 - 2;
                T x12 = 2 * C1_0;
                T x13 = x12 * x3;
                T x14 = 4 * C1_1 * pow(x3, 3);
                T x15 = 2 * x5;
                T x16 = pow(x2, -5.0 / 3.0);
                T x17 = x10 * x16;
                T x18 = (2.0 / 3.0) * x17;
                T x19 = f[0] * x15 - f[3] * x18;
                T x20 = 4 * x5;
                T x21 = (4.0 / 3.0) * x17;
                T x22 = f[0] * x20 - f[3] * x21;
                T x23 = C0_1 * x11;
                T x24 = f[1] * x15 + f[2] * x18;
                T x25 = f[1] * x20 + f[2] * x21;
                T x26 = f[1] * x18 + f[2] * x15;
                T x27 = f[1] * x21 + f[2] * x20;
                T x28 = -f[0] * x18 + f[3] * x15;
                T x29 = -f[0] * x21 + f[3] * x20;
                T x30 = 12 * C1_1 * x4;
                T x31 = x10 / pow(x2, 8.0 / 3.0);
                T x32 = x31 * x9;
                T x33 = (8.0 / 3.0) * x16;
                T x34 = -x0 * x33 + x15;
                T x35 = C0_1 * x19;
                T x36 = (16.0 / 3.0) * x16;
                T x37 = -x0 * x36 + x20;
                T x38 = (4.0 / 3.0) * x16;
                T x39 = f[0] * f[2];
                T x40 = x38 * x39;
                T x41 = f[1] * f[3];
                T x42 = x38 * x41;
                T x43 = f[2] * f[3];
                T x44 = x31 * x43;
                T x45 = x33 * x39;
                T x46 = x33 * x41;
                T x47 = C0_0 * (x40 - x42 - 10.0 / 9.0 * x44) - x12 * x43 + x23 * (-20.0 / 9.0 * x44 + x45 - x46) -
                        x30 * x43;
                T x48 = C0_1 * x24;
                T x49 = x31 * x8;
                T x50 = x1 * x33 + x15;
                T x51 = x1 * x36 + x20;
                T x52 = C0_1 * x26;
                T x53 = f[0] * f[1];
                T x54 = x38 * x53;
                T x55 = x38 * x43;
                T x56 = x31 * x41;
                T x57 = x33 * x53;
                T x58 = x33 * x43;
                T x59 = C0_0 * (x54 - x55 - 10.0 / 9.0 * x56) - x12 * x41 + x23 * (-20.0 / 9.0 * x56 + x57 - x58) -
                        x30 * x41;
                T x60 = x1 * x31;
                T x61 = C0_0 * (x18 + x38 * x7 + x38 * x8 + (10.0 / 9.0) * x60) + x1 * x12 + x1 * x30 - x13 - x14 +
                        x23 * (x21 + x33 * x7 + x33 * x8 + (20.0 / 9.0) * x60);
                T x62 = C0_1 * x28;
                T x63 = x31 * x39;
                T x64 = C0_0 * (-x54 + x55 - 10.0 / 9.0 * x63) - x12 * x39 + x23 * (-x57 + x58 - 20.0 / 9.0 * x63) -
                        x30 * x39;
                T x65 = x0 * x31;
                T x66 = C0_0 * (-x18 - x38 * x6 - x38 * x9 + (10.0 / 9.0) * x65) + x0 * x12 + x0 * x30 + x13 + x14 +
                        x23 * (-x21 - x33 * x6 - x33 * x9 + (20.0 / 9.0) * x65);
                T x67 = x31 * x7;
                T x68 = x31 * x53;
                T x69 = C0_0 * (-x40 + x42 - 10.0 / 9.0 * x68) - x12 * x53 + x23 * (-x45 + x46 - 20.0 / 9.0 * x68) -
                        x30 * x53;
                T x70 = x31 * x6;
                e += dx * (C0_0 * x11 + C0_1 * pow(x11, 2) + C1_0 * x4 + C1_1 * pow(x3, 4));
                lf[0] += dx * (grad_test[0] * (C0_0 * x19 + f[3] * x13 + f[3] * x14 + x22 * x23) +
                               grad_test[1] * (C0_0 * x24 - f[2] * x13 - f[2] * x14 + x23 * x25));
                lf[1] += dx * (grad_test[0] * (C0_0 * x26 - f[1] * x13 - f[1] * x14 + x23 * x27) +
                               grad_test[1] * (C0_0 * x28 + f[0] * x13 + f[0] * x14 + x23 * x29));
                bf[0] +=
                    dx * (grad_trial[0] * (grad_test[0] * (C0_0 * ((10.0 / 9.0) * x32 + x34) + x12 * x9 + x22 * x35 +
                                                           x23 * ((20.0 / 9.0) * x32 + x37) + x30 * x9) +
                                           grad_test[1] * (x25 * x35 + x47)) +
                          grad_trial[1] * (grad_test[0] * (x22 * x48 + x47) +
                                           grad_test[1] * (C0_0 * ((10.0 / 9.0) * x49 + x50) + x12 * x8 +
                                                           x23 * ((20.0 / 9.0) * x49 + x51) + x25 * x48 + x30 * x8)));
                bf[1] += dx * (grad_trial[0] * (grad_test[0] * (x22 * x52 + x59) + grad_test[1] * (x25 * x52 + x61)) +
                               grad_trial[1] * (grad_test[0] * (x22 * x62 + x66) + grad_test[1] * (x25 * x62 + x64)));
                bf[2] += dx * (grad_trial[0] * (grad_test[0] * (x27 * x35 + x59) + grad_test[1] * (x29 * x35 + x66)) +
                               grad_trial[1] * (grad_test[0] * (x27 * x48 + x61) + grad_test[1] * (x29 * x48 + x64)));
                bf[3] +=
                    dx * (grad_trial[0] * (grad_test[0] * (C0_0 * (x50 + (10.0 / 9.0) * x67) + x12 * x7 +
                                                           x23 * (x51 + (20.0 / 9.0) * x67) + x27 * x52 + x30 * x7) +
                                           grad_test[1] * (x29 * x52 + x69)) +
                          grad_trial[1] * (grad_test[0] * (x27 * x62 + x69) +
                                           grad_test[1] * (C0_0 * (x34 + (10.0 / 9.0) * x70) + x12 * x6 +
                                                           x23 * (x37 + (20.0 / 9.0) * x70) + x29 * x62 + x30 * x6)));
            }

            T C0_0{1.0};
            T C1_0{1.0};
            T C0_1{1.0};
            T C1_1{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_Yeoh_2_IMPL_hpp
