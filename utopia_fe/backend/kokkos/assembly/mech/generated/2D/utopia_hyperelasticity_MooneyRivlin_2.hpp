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
                }

                T C1{1.0};
                T C2{1.0};
            };

            MooneyRivlin(const Params &params) {
                C1 = params.C1;
                C2 = params.C2;
            }

            UTOPIA_FUNCTION void hessian(const T *UTOPIA_RESTRICT f,
                                         const T *grad_test,
                                         const T *grad_trial,
                                         const T dx,
                                         T *UTOPIA_RESTRICT bf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = 3 * C1;
                T x1 = grad_test[0] * grad_trial[0];
                T x2 = x0 * x1;
                T x3 = grad_test[1] * grad_trial[1];
                T x4 = x0 * x3;
                T x5 = C2 * x3;
                T x6 = pow(f[2], 3);
                T x7 = f[1] * x6;
                T x8 = C2 * x1;
                T x9 = pow(f[3], 3);
                T x10 = f[0] * x9;
                T x11 = grad_test[0] * grad_trial[1];
                T x12 = x0 * x11;
                T x13 = f[2] * x9;
                T x14 = grad_test[1] * grad_trial[0];
                T x15 = x0 * x14;
                T x16 = f[3] * x6;
                T x17 = pow(f[0], 2);
                T x18 = pow(f[3], 2);
                T x19 = x17 * x18;
                T x20 = pow(f[2], 2);
                T x21 = x18 * x2;
                T x22 = pow(f[1], 2);
                T x23 = x20 * x22;
                T x24 = x20 * x4;
                T x25 = 8 * C1;
                T x26 = x1 * x25;
                T x27 = x25 * x3;
                T x28 = f[2] * x11;
                T x29 = C1 * x17;
                T x30 = f[3] * x29;
                T x31 = x28 * x30;
                T x32 = C1 * x22;
                T x33 = f[3] * x32;
                T x34 = x28 * x33;
                T x35 = f[2] * x14;
                T x36 = x30 * x35;
                T x37 = x33 * x35;
                T x38 = f[0] * x18;
                T x39 = C2 * x38;
                T x40 = f[1] * f[2];
                T x41 = x18 * x40;
                T x42 = f[0] * f[3];
                T x43 = x20 * x42;
                T x44 = C2 * x11;
                T x45 = f[1] * x20;
                T x46 = f[3] * x45;
                T x47 = C2 * x14;
                T x48 = f[0] * x45;
                T x49 = C1 * x48;
                T x50 = x11 * x49;
                T x51 = C1 * x38;
                T x52 = f[1] * x51;
                T x53 = x11 * x52;
                T x54 = x14 * x49;
                T x55 = x14 * x52;
                T x56 = x40 * x42;
                T x57 = x17 * x24 + x21 * x22 - x26 * x56 - x27 * x56;
                T x58 = (1.0 / 4.0) * dx / (sqrt(-x40 + x42) * (x19 + x23 - 2 * x56));
                T x59 = C2 * x19;
                T x60 = x11 * x59;
                T x61 = C2 * x23;
                T x62 = x14 * x61;
                T x63 = pow(f[1], 3);
                T x64 = f[2] * x63;
                T x65 = pow(f[0], 3);
                T x66 = f[3] * x65;
                T x67 = x14 * x59;
                T x68 = 2 * C1;
                T x69 = x14 * x68;
                T x70 = x11 * x68;
                T x71 = x11 * x61;
                T x72 = x32 * x42;
                T x73 = x11 * x72;
                T x74 = C1 * x43;
                T x75 = x11 * x74;
                T x76 = x29 * x40;
                T x77 = x14 * x76;
                T x78 = C1 * x41;
                T x79 = x14 * x78;
                T x80 = x14 * x72;
                T x81 = x14 * x74;
                T x82 = x11 * x76;
                T x83 = x11 * x78;
                T x84 = 3 * x56;
                T x85 = x44 * x84;
                T x86 = x47 * x84;
                T x87 = f[2] * x3;
                T x88 = f[0] * x32;
                T x89 = f[1] * x30;
                T x90 = C1 * x46;
                T x91 = f[2] * f[3];
                T x92 = x17 * x5;
                T x93 = x22 * x8;
                T x94 = 4 * f[2] * x1;
                T x95 = 4 * x3;
                T x96 = -f[0] * x4 * x6 - f[1] * x2 * x9 + f[1] * x38 * x8 - f[2] * x4 * x65 - f[3] * x2 * x63 +
                        x1 * x89 + x1 * x90 - x48 * x5 + x51 * x87 - x51 * x94 + x87 * x88 - x88 * x94 - x89 * x95 -
                        x90 * x95 + x91 * x92 - x91 * x93;
                T x97 = f[0] * x63;
                T x98 = f[1] * x65;
                T x99 = x17 * x22;
                T x100 = f[1] * f[3] * x17;
                T x101 = C2 * f[0] * x22;
                bf[0] += x58 * (pow(f[2], 4) * x4 + pow(f[3], 4) * x2 - x10 * x8 - x12 * x13 - x12 * x16 - x13 * x15 -
                                x15 * x16 + x18 * x24 + x19 * x2 + x19 * x27 + x20 * x21 + x23 * x26 + x23 * x4 +
                                x28 * x39 + x31 + x34 + x35 * x39 + x36 + x37 + x41 * x8 - x43 * x5 - x44 * x46 -
                                x46 * x47 + x5 * x7 - 4 * x50 - 4 * x53 - 4 * x54 - 4 * x55 + x57);
                bf[1] += x58 * (-x10 * x12 + x10 * x69 - x12 * x66 - x15 * x64 - x15 * x7 + x60 - x62 + x64 * x70 +
                                x66 * x69 - 2 * x67 + x7 * x70 + 2 * x71 + x73 + x75 + x77 + x79 + 6 * x80 + 6 * x81 +
                                6 * x82 + 6 * x83 - x85 + x86 + x96);
                bf[2] += x58 * (-x10 * x15 + x10 * x70 - x12 * x64 - x12 * x7 - x15 * x66 - 2 * x60 + 2 * x62 +
                                x64 * x69 + x66 * x70 + x67 + x69 * x7 - x71 + 6 * x73 + 6 * x75 + 6 * x77 + 6 * x79 +
                                x80 + x81 + x82 + x83 + x85 - x86 + x96);
                bf[3] += x58 * (pow(f[0], 4) * x4 + pow(f[1], 4) * x2 + x100 * x44 + x100 * x47 - x101 * x28 -
                                x101 * x35 - x12 * x97 - x12 * x98 - x15 * x97 - x15 * x98 + x19 * x26 + x19 * x4 +
                                x2 * x23 + x2 * x99 + x23 * x27 - 4 * x31 - 4 * x34 - 4 * x36 - 4 * x37 + x4 * x99 +
                                x40 * x92 - x42 * x93 - x5 * x66 + x50 + x53 + x54 + x55 + x57 + x64 * x8);
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = C1 * grad_test[1];
                T x1 = C1 * grad_test[0];
                T x2 = pow(f[0], 2);
                T x3 = f[2] * x0;
                T x4 = pow(f[3], 2);
                T x5 = C2 * grad_test[0];
                T x6 = f[0] * x4;
                T x7 = C2 * grad_test[1];
                T x8 = pow(f[2], 2);
                T x9 = f[1] * x8;
                T x10 = f[1] * f[2];
                T x11 = f[0] * x10;
                T x12 = 4 * x1;
                T x13 = f[0] * f[3];
                T x14 = f[1] * x13;
                T x15 = 4 * x0;
                T x16 = f[2] * x13;
                T x17 = f[3] * x10;
                T x18 = f[3] * x1;
                T x19 = pow(f[1], 2);
                T x20 = (1.0 / 2.0) * dx / pow(-x10 + x13, 3.0 / 2.0);
                T x21 = f[1] * x1;
                T x22 = f[0] * x0;
                lf[0] +=
                    x20 * (pow(f[2], 3) * x0 - pow(f[3], 3) * x1 - x11 * x12 + x14 * x15 - x16 * x7 - x17 * x5 -
                           x18 * x19 + 3 * x18 * x2 - x18 * x8 - 3 * x19 * x3 + x2 * x3 + x3 * x4 + x5 * x6 + x7 * x9);
                lf[1] += x20 * (-pow(f[0], 3) * x0 + pow(f[1], 3) * x1 + f[2] * x19 * x5 + f[3] * x2 * x7 +
                                3 * x0 * x6 - 3 * x1 * x9 - x11 * x7 + x12 * x16 - x14 * x5 - x15 * x17 - x19 * x22 +
                                x2 * x21 + x21 * x4 - x22 * x8);
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = f[0] * f[3] - f[1] * f[2];
                T x1 = pow(f[0], 2) + pow(f[2], 2);
                T x2 = pow(f[1], 2) + pow(f[3], 2);
                T x3 = x1 + x2;
                e += dx * (C1 * (-2 + x3 / sqrt(x0)) +
                           C2 * (-2 + (-1.0 / 2.0 * pow(x1, 2) - 1.0 / 2.0 * pow(x2, 2) + (1.0 / 2.0) * pow(x3, 2) -
                                       pow(f[0] * f[1] + f[2] * f[3], 2)) /
                                          pow(x0, 3.0 / 2.0)));
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
                T x3 = pow(x2, -1.0 / 2.0);
                T x4 = pow(f[0], 2);
                T x5 = pow(f[2], 2);
                T x6 = x4 + x5;
                T x7 = pow(f[1], 2);
                T x8 = pow(f[3], 2);
                T x9 = x7 + x8;
                T x10 = x6 + x9;
                T x11 = pow(x2, -3.0 / 2.0);
                T x12 = f[0] * f[1];
                T x13 = f[2] * f[3];
                T x14 = pow(f[2], 3);
                T x15 = C1 * grad_test[1];
                T x16 = x14 * x15;
                T x17 = pow(f[3], 3);
                T x18 = C1 * grad_test[0];
                T x19 = x17 * x18;
                T x20 = f[2] * x15;
                T x21 = x20 * x8;
                T x22 = C2 * grad_test[0];
                T x23 = f[0] * x8;
                T x24 = x22 * x23;
                T x25 = C2 * grad_test[1];
                T x26 = f[1] * x5;
                T x27 = x25 * x26;
                T x28 = f[0] * x1;
                T x29 = 4 * x18;
                T x30 = f[1] * x0;
                T x31 = 4 * x15;
                T x32 = f[2] * x0;
                T x33 = f[3] * x1;
                T x34 = f[3] * x18;
                T x35 = 3 * x34;
                T x36 = 3 * x20;
                T x37 = (1.0 / 2.0) * dx * x11;
                T x38 = pow(f[1], 3);
                T x39 = x18 * x38;
                T x40 = pow(f[0], 3);
                T x41 = x15 * x40;
                T x42 = f[1] * x18;
                T x43 = f[3] * x25;
                T x44 = x4 * x43;
                T x45 = f[2] * x22;
                T x46 = x45 * x7;
                T x47 = f[0] * x15;
                T x48 = 3 * x15;
                T x49 = 3 * x18;
                T x50 = grad_trial[0] * x49;
                T x51 = grad_trial[1] * x48;
                T x52 = grad_trial[1] * x14;
                T x53 = grad_trial[0] * x17;
                T x54 = f[2] * grad_trial[1];
                T x55 = 3 * x19;
                T x56 = f[3] * grad_trial[0];
                T x57 = 3 * x16;
                T x58 = x4 * x8;
                T x59 = x50 * x8;
                T x60 = x5 * x7;
                T x61 = x5 * x51;
                T x62 = 8 * x60;
                T x63 = grad_trial[0] * x18;
                T x64 = 8 * x58;
                T x65 = grad_trial[1] * x15;
                T x66 = x13 * x4;
                T x67 = grad_trial[1] * x18;
                T x68 = x13 * x7;
                T x69 = grad_trial[0] * x15;
                T x70 = grad_trial[0] * x25;
                T x71 = f[2] * x70;
                T x72 = x1 * x8;
                T x73 = grad_trial[0] * x22;
                T x74 = x0 * x5;
                T x75 = grad_trial[1] * x25;
                T x76 = grad_trial[1] * x22;
                T x77 = f[3] * x76;
                T x78 = x12 * x5;
                T x79 = grad_trial[1] * x29;
                T x80 = x12 * x8;
                T x81 = grad_trial[0] * x31;
                T x82 = x0 * x1;
                T x83 = 8 * x82;
                T x84 = x4 * x61 + x59 * x7 - x63 * x83 - x65 * x83;
                T x85 = (1.0 / 4.0) * dx * x3 / (x58 + x60 - 2 * x82);
                T x86 = x58 * x76;
                T x87 = x60 * x70;
                T x88 = f[0] * grad_trial[1];
                T x89 = f[1] * grad_trial[0];
                T x90 = grad_trial[0] * x38;
                T x91 = grad_trial[1] * x40;
                T x92 = x58 * x70;
                T x93 = x47 * x53;
                T x94 = x42 * x52;
                T x95 = x60 * x76;
                T x96 = x0 * x7;
                T x97 = x67 * x96;
                T x98 = x67 * x74;
                T x99 = x1 * x4;
                T x100 = x69 * x99;
                T x101 = x69 * x72;
                T x102 = x69 * x96;
                T x103 = x69 * x74;
                T x104 = x67 * x99;
                T x105 = x67 * x72;
                T x106 = 3 * x82;
                T x107 = x106 * x76;
                T x108 = x106 * x70;
                T x109 = 3 * x54;
                T x110 = 3 * x56;
                T x111 = f[0] * x7;
                T x112 = f[2] * grad_trial[0] * x29;
                T x113 = f[1] * x4;
                T x114 = f[3] * grad_trial[1] * x31;
                T x115 = grad_trial[0] * x26 * x34 - x109 * x41 - x110 * x39 - x111 * x112 - x112 * x23 - x113 * x114 -
                         x114 * x26 + x20 * x7 * x88 + x21 * x88 + x34 * x4 * x89 - x55 * x89 - x57 * x88 + x66 * x75 -
                         x68 * x73 + x73 * x80 - x75 * x78;
                T x116 = x4 * x7;
                e += dx * (C1 * (x10 * x3 - 2) + C2 * (x11 * ((1.0 / 2.0) * pow(x10, 2) - 1.0 / 2.0 * pow(x6, 2) -
                                                              1.0 / 2.0 * pow(x9, 2) - pow(x12 + x13, 2)) -
                                                       2));
                lf[0] += x37 * (x16 - x19 + x20 * x4 + x21 - x22 * x33 + x24 - x25 * x32 + x27 - x28 * x29 + x30 * x31 -
                                x34 * x5 - x34 * x7 + x35 * x4 - x36 * x7);
                lf[1] += x37 * (-x22 * x30 + x23 * x48 - x25 * x28 - x26 * x49 + x29 * x32 - x31 * x33 + x39 +
                                x4 * x42 - x41 + x42 * x8 + x44 + x46 - x47 * x5 - x47 * x7);
                bf[0] += x85 *
                         (-f[0] * x22 * x53 + f[1] * x25 * x52 + pow(f[2], 4) * x51 + pow(f[3], 4) * x50 + x23 * x71 +
                          x24 * x54 - x26 * x77 - x27 * x56 - x35 * x52 - x36 * x53 + x5 * x59 + x50 * x58 + x51 * x60 -
                          x54 * x55 - x56 * x57 + x61 * x8 + x62 * x63 + x64 * x65 + x66 * x67 + x66 * x69 + x67 * x68 +
                          x68 * x69 + x72 * x73 - x74 * x75 - x78 * x79 - x78 * x81 - x79 * x80 - x80 * x81 + x84);
                bf[1] += x85 * (x100 + x101 + 6 * x102 + 6 * x103 + 6 * x104 + 6 * x105 - x107 + x108 + x115 -
                                x35 * x91 - x36 * x90 + 2 * x39 * x54 + 2 * x41 * x56 - x55 * x88 - x57 * x89 + x86 -
                                x87 - 2 * x92 + 2 * x93 + 2 * x94 + 2 * x95 + x97 + x98);
                bf[2] += x85 * (6 * x100 + 6 * x101 + x102 + x103 + x104 + x105 + x107 - x108 - x109 * x39 -
                                x110 * x41 + x115 + 2 * x16 * x89 + 2 * x19 * x88 + 2 * x20 * x90 + 2 * x34 * x91 -
                                2 * x86 + 2 * x87 + x92 - 3 * x93 - 3 * x94 - x95 + 6 * x97 + 6 * x98);
                bf[3] += x85 * (pow(f[0], 4) * x51 + pow(f[1], 4) * x50 - x111 * x71 + x113 * x77 + x116 * x50 +
                                x116 * x51 - 3 * x39 * x88 - 3 * x41 * x89 - 3 * x42 * x91 - x43 * x91 + x44 * x89 +
                                x45 * x90 - x46 * x88 - 3 * x47 * x90 + x50 * x60 + x51 * x58 + x62 * x65 + x63 * x64 -
                                x66 * x79 - x66 * x81 + x67 * x78 + x67 * x80 - x68 * x79 - x68 * x81 + x69 * x78 +
                                x69 * x80 - x73 * x96 + x75 * x99 + x84);
            }

            T C1{1.0};
            T C2{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_MooneyRivlin_2_IMPL_hpp
