// #ifndef UTOPIA_TPL_HYPERELASTICITY_NeoHookeanSiguenza_2_IMPL_hpp
// #define UTOPIA_TPL_HYPERELASTICITY_NeoHookeanSiguenza_2_IMPL_hpp

// #include "utopia_Input.hpp"

// #include "utopia_hyperelasticity_NeoHookeanSiguenza.hpp"

// namespace utopia {
//     namespace kernels {

//         /**
//          * Specialization of NeoHookeanSiguenza for dimension 2
//          */
//         template <typename T>
//         class NeoHookeanSiguenza<T, 2> {
//         public:
//             static constexpr int Dim = 2;

//             UTOPIA_FUNCTION static constexpr const char *class_name() { return "NeoHookeanSiguenza_2"; }

//             class Params : public Configurable {
//             public:
//                 void read(Input &in) override {
//                     in.get("G", G);
//                     in.get("K", K);
//                 }

//                 T G{1.0};
//                 T K{1.0};
//             };

//             NeoHookeanSiguenza(const Params &params) {
//                 G = params.G;
//                 K = params.K;
//             }

//             UTOPIA_FUNCTION void hessian(const T *UTOPIA_RESTRICT f,
//                                          const T *grad_test,
//                                          const T *grad_trial,
//                                          const T dx,
//                                          T *UTOPIA_RESTRICT bf) const {
//                 using namespace utopia::device;
//                 // Automatically generated
//                 T x0 = f[0] * f[3];
//                 T x1 = f[1] * f[2];
//                 T x2 = x0 - x1;
//                 T x3 = K / pow(x2, 2);
//                 T x4 = f[3] * x3;
//                 T x5 = f[2] * x4;
//                 T x6 = log(x2);
//                 T x7 = pow(x2, -1.6666666666666665);
//                 T x8 = 1.3333333333333333 * x7;
//                 T x9 = f[0] * x8;
//                 T x10 = f[2] * x9;
//                 T x11 = f[3] * x8;
//                 T x12 = f[1] * x11;
//                 T x13 = pow(f[0], 2);
//                 T x14 = pow(f[1], 2);
//                 T x15 = pow(f[2], 2);
//                 T x16 = pow(f[3], 2);
//                 T x17 = x13 + x14 + x15 + x16;
//                 T x18 = pow(x2, -2.6666666666666665);
//                 T x19 = 1.1111111111111109 * x17 * x18;
//                 T x20 = f[3] * x19;
//                 T x21 = (1.0 / 2.0) * G;
//                 T x22 = x21 * (-f[2] * x20 + x10 - x12) + x5 * x6 - x5;
//                 T x23 = x16 * x3;
//                 T x24 = 2 * pow(x2, -0.66666666666666663);
//                 T x25 = 2.6666666666666665 * x7;
//                 T x26 = -x0 * x25 + x24;
//                 T x27 = x15 * x3;
//                 T x28 = x1 * x25 + x24;
//                 T x29 = f[1] * x4;
//                 T x30 = f[1] * x9;
//                 T x31 = f[2] * x11;
//                 T x32 = grad_trial[0] * (x21 * (-f[1] * x20 + x30 - x31) + x29 * x6 - x29);
//                 T x33 = x0 * x3;
//                 T x34 = K * x6 / x2;
//                 T x35 = 0.66666666666666663 * x17 * x7;
//                 T x36 = x21 * (1.1111111111111109 * f[0] * f[3] * x17 * x18 - x13 * x8 - x16 * x8 - x35) - x33 * x6 +
//                         x33 + x34;
//                 T x37 = f[0] * x3;
//                 T x38 = f[2] * x37;
//                 T x39 = f[0] * x19;
//                 T x40 = grad_trial[1] * (x21 * (-f[2] * x39 - x30 + x31) + x38 * x6 - x38);
//                 T x41 = x1 * x3;
//                 T x42 = x21 * (x1 * x19 + x14 * x8 + x15 * x8 + x35) - x34 - x41 * x6 + x41;
//                 T x43 = f[1] * x37;
//                 T x44 = x21 * (-f[1] * x39 - x10 + x12) + x43 * x6 - x43;
//                 T x45 = x14 * x3;
//                 T x46 = x13 * x3;
//                 bf[0] +=
//                     dx *
//                     (grad_test[0] * (grad_trial[0] * (x21 * (x16 * x19 + x26) - x23 * x6 + x23) + grad_trial[1] *
//                     x22) +
//                      grad_test[1] * (grad_trial[0] * x22 + grad_trial[1] * (x21 * (x15 * x19 + x28) - x27 * x6 +
//                      x27)));
//                 bf[1] += dx * (grad_test[0] * (grad_trial[1] * x36 + x32) + grad_test[1] * (grad_trial[0] * x42 +
//                 x40)); bf[2] += dx * (grad_test[0] * (grad_trial[1] * x42 + x32) + grad_test[1] * (grad_trial[0] *
//                 x36 + x40)); bf[3] +=
//                     dx *
//                     (grad_test[0] * (grad_trial[0] * (x21 * (x14 * x19 + x28) - x45 * x6 + x45) + grad_trial[1] *
//                     x44) +
//                      grad_test[1] * (grad_trial[0] * x44 + grad_trial[1] * (x21 * (x13 * x19 + x26) - x46 * x6 +
//                      x46)));
//             }

//             UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
//                                           const T *UTOPIA_RESTRICT grad_test,
//                                           const T dx,
//                                           T *UTOPIA_RESTRICT lf) const {
//                 using namespace utopia::device;
//                 // Automatically generated
//                 T x0 = f[0] * f[3] - f[1] * f[2];
//                 T x1 = K * log(x0) / x0;
//                 T x2 = pow(x0, -0.66666666666666663);
//                 T x3 = 2 * x2;
//                 T x4 = 0.66666666666666663 * pow(x0, -1.6666666666666665) *
//                        (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2));
//                 T x5 = (1.0 / 2.0) * G;
//                 lf[0] += dx * (grad_test[0] * (f[3] * x1 + x5 * (f[0] * x3 - f[3] * x4)) +
//                                grad_test[1] * (-f[2] * x1 + x5 * (f[1] * x3 + f[2] * x4)));
//                 lf[1] += dx * (grad_test[0] * (-f[1] * x1 + x5 * (f[1] * x4 + f[2] * x3)) +
//                                grad_test[1] * (f[0] * x1 + x5 * (-f[0] * x4 + 2 * f[3] * x2)));
//             }

//             UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
//                 using namespace utopia::device;
//                 // Automatically generated
//                 T x0 = f[0] * f[3] - f[1] * f[2];
//                 e += dx *
//                      ((1.0 / 2.0) * G *
//                           (pow(x0, -0.66666666666666663) * (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3],
//                           2)) -
//                            2) +
//                       (1.0 / 2.0) * K * pow(log(x0), 2));
//             }

//             UTOPIA_FUNCTION void eval(const T *UTOPIA_RESTRICT f,
//                                       const T *grad_test,
//                                       const T *grad_trial,
//                                       const T dx,
//                                       T &e,
//                                       T *UTOPIA_RESTRICT lf,
//                                       T *UTOPIA_RESTRICT bf) const {
//                 using namespace utopia::device;
//                 // Automatically generated
//                 T x0 = f[0] * f[3];
//                 T x1 = f[1] * f[2];
//                 T x2 = x0 - x1;
//                 T x3 = log(x2);
//                 T x4 = pow(x2, -0.66666666666666663);
//                 T x5 = pow(f[0], 2);
//                 T x6 = pow(f[1], 2);
//                 T x7 = pow(f[2], 2);
//                 T x8 = pow(f[3], 2);
//                 T x9 = x5 + x6 + x7 + x8;
//                 T x10 = (1.0 / 2.0) * G;
//                 T x11 = K * x3 / x2;
//                 T x12 = 2 * x4;
//                 T x13 = pow(x2, -1.6666666666666665);
//                 T x14 = 0.66666666666666663 * x13 * x9;
//                 T x15 = K / pow(x2, 2);
//                 T x16 = f[3] * x15;
//                 T x17 = f[2] * x16;
//                 T x18 = 1.3333333333333333 * x13;
//                 T x19 = f[0] * x18;
//                 T x20 = f[2] * x19;
//                 T x21 = f[3] * x18;
//                 T x22 = f[1] * x21;
//                 T x23 = pow(x2, -2.6666666666666665);
//                 T x24 = 1.1111111111111109 * x23 * x9;
//                 T x25 = f[3] * x24;
//                 T x26 = x10 * (-f[2] * x25 + x20 - x22) + x17 * x3 - x17;
//                 T x27 = x15 * x8;
//                 T x28 = 2.6666666666666665 * x13;
//                 T x29 = -x0 * x28 + x12;
//                 T x30 = x15 * x7;
//                 T x31 = x1 * x28 + x12;
//                 T x32 = f[1] * x16;
//                 T x33 = f[1] * x19;
//                 T x34 = f[2] * x21;
//                 T x35 = grad_trial[0] * (x10 * (-f[1] * x25 + x33 - x34) + x3 * x32 - x32);
//                 T x36 = x0 * x15;
//                 T x37 = x10 * (1.1111111111111109 * f[0] * f[3] * x23 * x9 - x14 - x18 * x5 - x18 * x8) + x11 -
//                         x3 * x36 + x36;
//                 T x38 = f[0] * x15;
//                 T x39 = f[2] * x38;
//                 T x40 = f[0] * x24;
//                 T x41 = grad_trial[1] * (x10 * (-f[2] * x40 - x33 + x34) + x3 * x39 - x39);
//                 T x42 = x1 * x15;
//                 T x43 = x10 * (x1 * x24 + x14 + x18 * x6 + x18 * x7) - x11 - x3 * x42 + x42;
//                 T x44 = f[1] * x38;
//                 T x45 = x10 * (-f[1] * x40 - x20 + x22) + x3 * x44 - x44;
//                 T x46 = x15 * x6;
//                 T x47 = x15 * x5;
//                 e += dx * ((1.0 / 2.0) * K * pow(x3, 2) + x10 * (x4 * x9 - 2));
//                 lf[0] += dx * (grad_test[0] * (f[3] * x11 + x10 * (f[0] * x12 - f[3] * x14)) +
//                                grad_test[1] * (-f[2] * x11 + x10 * (f[1] * x12 + f[2] * x14)));
//                 lf[1] += dx * (grad_test[0] * (-f[1] * x11 + x10 * (f[1] * x14 + f[2] * x12)) +
//                                grad_test[1] * (f[0] * x11 + x10 * (-f[0] * x14 + 2 * f[3] * x4)));
//                 bf[0] +=
//                     dx *
//                     (grad_test[0] * (grad_trial[0] * (x10 * (x24 * x8 + x29) - x27 * x3 + x27) + grad_trial[1] * x26)
//                     +
//                      grad_test[1] * (grad_trial[0] * x26 + grad_trial[1] * (x10 * (x24 * x7 + x31) - x3 * x30 +
//                      x30)));
//                 bf[1] += dx * (grad_test[0] * (grad_trial[1] * x37 + x35) + grad_test[1] * (grad_trial[0] * x43 +
//                 x41)); bf[2] += dx * (grad_test[0] * (grad_trial[1] * x43 + x35) + grad_test[1] * (grad_trial[0] *
//                 x37 + x41)); bf[3] +=
//                     dx *
//                     (grad_test[0] * (grad_trial[0] * (x10 * (x24 * x6 + x31) - x3 * x46 + x46) + grad_trial[1] * x45)
//                     +
//                      grad_test[1] * (grad_trial[0] * x45 + grad_trial[1] * (x10 * (x24 * x5 + x29) - x3 * x47 +
//                      x47)));
//             }

//             UTOPIA_FUNCTION void apply(const T *UTOPIA_RESTRICT f,
//                                        const T *grad_test,
//                                        const T *disp_grad,
//                                        const T dx,
//                                        T *UTOPIA_RESTRICT res) const {
//                 using namespace utopia::device;
//                 // Automatically generated
//                 T x0 = f[0] * f[3];
//                 T x1 = f[1] * f[2];
//                 T x2 = x0 - x1;
//                 T x3 = K / pow(x2, 2);
//                 T x4 = f[3] * x3;
//                 T x5 = f[2] * x4;
//                 T x6 = log(x2);
//                 T x7 = pow(x2, -1.6666666666666665);
//                 T x8 = 1.3333333333333333 * x7;
//                 T x9 = f[0] * x8;
//                 T x10 = f[2] * x9;
//                 T x11 = f[3] * x8;
//                 T x12 = f[1] * x11;
//                 T x13 = pow(f[0], 2);
//                 T x14 = pow(f[1], 2);
//                 T x15 = pow(f[2], 2);
//                 T x16 = pow(f[3], 2);
//                 T x17 = x13 + x14 + x15 + x16;
//                 T x18 = pow(x2, -2.6666666666666665);
//                 T x19 = 1.1111111111111109 * x17 * x18;
//                 T x20 = f[3] * x19;
//                 T x21 = (1.0 / 2.0) * G;
//                 T x22 = x21 * (-f[2] * x20 + x10 - x12) + x5 * x6 - x5;
//                 T x23 = f[1] * x4;
//                 T x24 = f[1] * x9;
//                 T x25 = f[2] * x11;
//                 T x26 = x21 * (-f[1] * x20 + x24 - x25) + x23 * x6 - x23;
//                 T x27 = x16 * x3;
//                 T x28 = 2 * pow(x2, -0.66666666666666663);
//                 T x29 = 2.6666666666666665 * x7;
//                 T x30 = -x0 * x29 + x28;
//                 T x31 = x0 * x3;
//                 T x32 = K * x6 / x2;
//                 T x33 = 0.66666666666666663 * x17 * x7;
//                 T x34 = x21 * (1.1111111111111109 * f[0] * f[3] * x17 * x18 - x13 * x8 - x16 * x8 - x33) - x31 * x6 +
//                         x31 + x32;
//                 T x35 = f[0] * x3;
//                 T x36 = f[2] * x35;
//                 T x37 = f[0] * x19;
//                 T x38 = x21 * (-f[2] * x37 - x24 + x25) + x36 * x6 - x36;
//                 T x39 = x15 * x3;
//                 T x40 = x1 * x29 + x28;
//                 T x41 = x1 * x3;
//                 T x42 = x21 * (x1 * x19 + x14 * x8 + x15 * x8 + x33) - x32 - x41 * x6 + x41;
//                 T x43 = f[1] * x35;
//                 T x44 = x21 * (-f[1] * x37 - x10 + x12) + x43 * x6 - x43;
//                 T x45 = x13 * x3;
//                 T x46 = x14 * x3;
//                 res[0] += dx * (grad_test[0] * (disp_grad[0] * (x21 * (x16 * x19 + x30) - x27 * x6 + x27) +
//                                                 disp_grad[1] * x22 + disp_grad[2] * x26 + disp_grad[3] * x34) +
//                                 grad_test[1] *
//                                     (disp_grad[0] * x22 + disp_grad[1] * (x21 * (x15 * x19 + x40) - x39 * x6 + x39) +
//                                      disp_grad[2] * x42 + disp_grad[3] * x38));
//                 res[1] += dx * (grad_test[0] *
//                                     (disp_grad[0] * x26 + disp_grad[1] * x42 +
//                                      disp_grad[2] * (x21 * (x14 * x19 + x40) - x46 * x6 + x46) + disp_grad[3] * x44)
//                                      +
//                                 grad_test[1] * (disp_grad[0] * x34 + disp_grad[1] * x38 + disp_grad[2] * x44 +
//                                                 disp_grad[3] * (x21 * (x13 * x19 + x30) - x45 * x6 + x45)));
//             }

//             T G{1.0};
//             T K{1.0};
//         };
//     }  // namespace kernels
// }  // namespace utopia

// #endif  // UTOPIA_TPL_HYPERELASTICITY_NeoHookeanSiguenza_2_IMPL_hpp
