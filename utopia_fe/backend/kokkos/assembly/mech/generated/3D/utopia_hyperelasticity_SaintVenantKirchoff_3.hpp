#ifndef UTOPIA_TPL_HYPERELASTICITY_SaintVenantKirchoff_3_IMPL_hpp
#define UTOPIA_TPL_HYPERELASTICITY_SaintVenantKirchoff_3_IMPL_hpp

#include "utopia_Input.hpp"

#include "utopia_hyperelasticity_SaintVenantKirchoff.hpp"

namespace utopia {
    namespace kernels {

        /**
         * Specialization of SaintVenantKirchoff for dimension 3
         */
        template <typename T>
        class SaintVenantKirchoff<T, 3> {
        public:
            static constexpr int Dim = 3;

            UTOPIA_FUNCTION static constexpr const char *class_name() { return "SaintVenantKirchoff_3"; }

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
                T x1 = f[3] * f[4];
                T x2 = f[6] * f[7];
                T x3 = lambda * x0 + mu * (2 * x0 + x1 + x2);
                T x4 = f[0] * f[2];
                T x5 = f[3] * f[5];
                T x6 = f[6] * f[8];
                T x7 = lambda * x4 + mu * (2 * x4 + x5 + x6);
                T x8 = pow(f[0], 2);
                T x9 = pow(f[1], 2);
                T x10 = pow(f[3], 2);
                T x11 = x10 + x9;
                T x12 = pow(f[6], 2);
                T x13 = pow(f[2], 2);
                T x14 = x13 - 1;
                T x15 = x12 + x14;
                T x16 = pow(f[4], 2);
                T x17 = pow(f[5], 2);
                T x18 = pow(f[7], 2);
                T x19 = pow(f[8], 2);
                T x20 = lambda * ((1.0 / 2.0) * x10 + (1.0 / 2.0) * x12 + (1.0 / 2.0) * x13 + (1.0 / 2.0) * x16 +
                                  (1.0 / 2.0) * x17 + (1.0 / 2.0) * x18 + (1.0 / 2.0) * x19 + (1.0 / 2.0) * x8 +
                                  (1.0 / 2.0) * x9 - 3.0 / 2.0);
                T x21 = f[1] * f[2];
                T x22 = f[4] * f[5];
                T x23 = f[7] * f[8];
                T x24 = lambda * x21 + mu * (2 * x21 + x22 + x23);
                T x25 = x16 + x8;
                T x26 = x17 - 1;
                T x27 = x19 + x8;
                T x28 = f[0] * mu;
                T x29 = f[3] * lambda;
                T x30 = f[1] * x29 + f[4] * x28;
                T x31 = f[2] * x29 + f[5] * x28;
                T x32 = f[0] * f[3];
                T x33 = f[1] * f[4];
                T x34 = f[2] * f[5];
                T x35 = grad_test[0] * (lambda * x32 + mu * (2 * x32 + x33 + x34));
                T x36 = f[4] * lambda;
                T x37 = f[1] * mu;
                T x38 = f[0] * x36 + f[3] * x37;
                T x39 = f[2] * x36 + f[5] * x37;
                T x40 = grad_test[1] * (lambda * x33 + mu * (x32 + 2 * x33 + x34));
                T x41 = f[5] * lambda;
                T x42 = f[2] * mu;
                T x43 = f[0] * x41 + f[3] * x42;
                T x44 = f[1] * x41 + f[4] * x42;
                T x45 = grad_test[2] * (lambda * x34 + mu * (x32 + x33 + 2 * x34));
                T x46 = f[6] * lambda;
                T x47 = f[1] * x46 + f[7] * x28;
                T x48 = f[2] * x46 + f[8] * x28;
                T x49 = f[0] * f[6];
                T x50 = f[1] * f[7];
                T x51 = f[2] * f[8];
                T x52 = grad_test[0] * (lambda * x49 + mu * (2 * x49 + x50 + x51));
                T x53 = f[7] * lambda;
                T x54 = f[0] * x53 + f[6] * x37;
                T x55 = f[2] * x53 + f[8] * x37;
                T x56 = grad_test[1] * (lambda * x50 + mu * (x49 + 2 * x50 + x51));
                T x57 = f[8] * lambda;
                T x58 = f[0] * x57 + f[6] * x42;
                T x59 = f[1] * x57 + f[7] * x42;
                T x60 = grad_test[2] * (lambda * x51 + mu * (x49 + x50 + 2 * x51));
                T x61 = lambda * x1 + mu * (x0 + 2 * x1 + x2);
                T x62 = lambda * x5 + mu * (x4 + 2 * x5 + x6);
                T x63 = lambda * x22 + mu * (x21 + 2 * x22 + x23);
                T x64 = x16 + x19;
                T x65 = f[3] * mu;
                T x66 = f[6] * x36 + f[7] * x65;
                T x67 = f[6] * x41 + f[8] * x65;
                T x68 = f[3] * f[6];
                T x69 = f[4] * f[7];
                T x70 = f[5] * f[8];
                T x71 = grad_test[0] * (lambda * x68 + mu * (2 * x68 + x69 + x70));
                T x72 = f[4] * mu;
                T x73 = f[6] * x72 + f[7] * x29;
                T x74 = f[7] * x41 + f[8] * x72;
                T x75 = grad_test[1] * (lambda * x69 + mu * (x68 + 2 * x69 + x70));
                T x76 = f[5] * mu;
                T x77 = f[6] * x76 + f[8] * x29;
                T x78 = f[7] * x76 + f[8] * x36;
                T x79 = grad_test[2] * (lambda * x70 + mu * (x68 + x69 + 2 * x70));
                T x80 = lambda * x2 + mu * (x0 + x1 + 2 * x2);
                T x81 = lambda * x6 + mu * (x4 + x5 + 2 * x6);
                T x82 = lambda * x23 + mu * (x21 + x22 + 2 * x23);
                bf[0] += dx * (grad_trial[0] * (grad_test[0] * (lambda * x8 + mu * (x11 + x15 + 3 * x8) + x20) +
                                                grad_test[1] * x3 + grad_test[2] * x7) +
                               grad_trial[1] * (grad_test[0] * x3 +
                                                grad_test[1] * (lambda * x9 + mu * (x14 + x18 + x25 + 3 * x9) + x20) +
                                                grad_test[2] * x24) +
                               grad_trial[2] * (grad_test[0] * x7 + grad_test[1] * x24 +
                                                grad_test[2] * (lambda * x13 + mu * (3 * x13 + x26 + x27 + x9) + x20)));
                bf[1] += dx * (grad_trial[0] * (grad_test[1] * x30 + grad_test[2] * x31 + x35) +
                               grad_trial[1] * (grad_test[0] * x38 + grad_test[2] * x39 + x40) +
                               grad_trial[2] * (grad_test[0] * x43 + grad_test[1] * x44 + x45));
                bf[2] += dx * (grad_trial[0] * (grad_test[1] * x47 + grad_test[2] * x48 + x52) +
                               grad_trial[1] * (grad_test[0] * x54 + grad_test[2] * x55 + x56) +
                               grad_trial[2] * (grad_test[0] * x58 + grad_test[1] * x59 + x60));
                bf[3] += dx * (grad_trial[0] * (grad_test[1] * x38 + grad_test[2] * x43 + x35) +
                               grad_trial[1] * (grad_test[0] * x30 + grad_test[2] * x44 + x40) +
                               grad_trial[2] * (grad_test[0] * x31 + grad_test[1] * x39 + x45));
                bf[4] +=
                    dx * (grad_trial[0] * (grad_test[0] * (lambda * x10 + mu * (3 * x10 + x12 + x25 + x26) + x20) +
                                           grad_test[1] * x61 + grad_test[2] * x62) +
                          grad_trial[1] * (grad_test[0] * x61 +
                                           grad_test[1] * (lambda * x16 + mu * (x11 + 3 * x16 + x18 + x26) + x20) +
                                           grad_test[2] * x63) +
                          grad_trial[2] * (grad_test[0] * x62 + grad_test[1] * x63 +
                                           grad_test[2] * (lambda * x17 + mu * (x10 + x14 + 3 * x17 + x64) + x20)));
                bf[5] += dx * (grad_trial[0] * (grad_test[1] * x66 + grad_test[2] * x67 + x71) +
                               grad_trial[1] * (grad_test[0] * x73 + grad_test[2] * x74 + x75) +
                               grad_trial[2] * (grad_test[0] * x77 + grad_test[1] * x78 + x79));
                bf[6] += dx * (grad_trial[0] * (grad_test[1] * x54 + grad_test[2] * x58 + x52) +
                               grad_trial[1] * (grad_test[0] * x47 + grad_test[2] * x59 + x56) +
                               grad_trial[2] * (grad_test[0] * x48 + grad_test[1] * x55 + x60));
                bf[7] += dx * (grad_trial[0] * (grad_test[1] * x73 + grad_test[2] * x77 + x71) +
                               grad_trial[1] * (grad_test[0] * x66 + grad_test[2] * x78 + x75) +
                               grad_trial[2] * (grad_test[0] * x67 + grad_test[1] * x74 + x79));
                bf[8] +=
                    dx * (grad_trial[0] * (grad_test[0] * (lambda * x12 + mu * (x10 + 3 * x12 + x18 + x27 - 1) + x20) +
                                           grad_test[1] * x80 + grad_test[2] * x81) +
                          grad_trial[1] * (grad_test[0] * x80 +
                                           grad_test[1] * (lambda * x18 + mu * (x12 + 3 * x18 + x64 + x9 - 1) + x20) +
                                           grad_test[2] * x82) +
                          grad_trial[2] * (grad_test[0] * x81 + grad_test[1] * x82 +
                                           grad_test[2] * (lambda * x19 + mu * (x15 + x17 + x18 + 3 * x19) + x20)));
            }

            UTOPIA_FUNCTION void gradient(const T *UTOPIA_RESTRICT f,
                                          const T *UTOPIA_RESTRICT grad_test,
                                          const T dx,
                                          T *UTOPIA_RESTRICT lf) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = (1.0 / 2.0) * pow(f[0], 2) + (1.0 / 2.0) * pow(f[3], 2) + (1.0 / 2.0) * pow(f[6], 2);
                T x1 = (1.0 / 2.0) * pow(f[1], 2) + (1.0 / 2.0) * pow(f[4], 2) + (1.0 / 2.0) * pow(f[7], 2);
                T x2 = (1.0 / 2.0) * pow(f[2], 2) + (1.0 / 2.0) * pow(f[5], 2) + (1.0 / 2.0) * pow(f[8], 2);
                T x3 = lambda * (x0 + x1 + x2 - 3.0 / 2.0);
                T x4 = (1.0 / 2.0) * f[0];
                T x5 = (1.0 / 2.0) * f[3];
                T x6 = (1.0 / 2.0) * f[6];
                T x7 = f[1] * x4 + f[4] * x5 + f[7] * x6;
                T x8 = 2 * f[1];
                T x9 = f[2] * x4 + f[5] * x5 + f[8] * x6;
                T x10 = 2 * f[2];
                T x11 = x0 - 1.0 / 2.0;
                T x12 = 2 * f[0];
                T x13 = (1.0 / 2.0) * f[1] * f[2] + (1.0 / 2.0) * f[4] * f[5] + (1.0 / 2.0) * f[7] * f[8];
                T x14 = x1 - 1.0 / 2.0;
                T x15 = x2 - 1.0 / 2.0;
                T x16 = 2 * f[4];
                T x17 = 2 * f[5];
                T x18 = 2 * f[3];
                T x19 = 2 * f[7];
                T x20 = 2 * f[8];
                T x21 = 2 * f[6];
                lf[0] += dx * (grad_test[0] * (f[0] * x3 + mu * (x10 * x9 + x11 * x12 + x7 * x8)) +
                               grad_test[1] * (f[1] * x3 + mu * (x10 * x13 + x12 * x7 + x14 * x8)) +
                               grad_test[2] * (f[2] * x3 + mu * (x10 * x15 + x12 * x9 + x13 * x8)));
                lf[1] += dx * (grad_test[0] * (f[3] * x3 + mu * (x11 * x18 + x16 * x7 + x17 * x9)) +
                               grad_test[1] * (f[4] * x3 + mu * (x13 * x17 + x14 * x16 + x18 * x7)) +
                               grad_test[2] * (f[5] * x3 + mu * (x13 * x16 + x15 * x17 + x18 * x9)));
                lf[2] += dx * (grad_test[0] * (f[6] * x3 + mu * (x11 * x21 + x19 * x7 + x20 * x9)) +
                               grad_test[1] * (f[7] * x3 + mu * (x13 * x20 + x14 * x19 + x21 * x7)) +
                               grad_test[2] * (f[8] * x3 + mu * (x13 * x19 + x15 * x20 + x21 * x9)));
            }

            UTOPIA_FUNCTION void value(const T *UTOPIA_RESTRICT f, const T dx, T &e) const {
                using namespace utopia::device;
                // Automatically generated
                T x0 = (1.0 / 2.0) * pow(f[0], 2) + (1.0 / 2.0) * pow(f[3], 2) + (1.0 / 2.0) * pow(f[6], 2);
                T x1 = (1.0 / 2.0) * pow(f[1], 2) + (1.0 / 2.0) * pow(f[4], 2) + (1.0 / 2.0) * pow(f[7], 2);
                T x2 = (1.0 / 2.0) * pow(f[2], 2) + (1.0 / 2.0) * pow(f[5], 2) + (1.0 / 2.0) * pow(f[8], 2);
                T x3 = (1.0 / 2.0) * f[0];
                T x4 = (1.0 / 2.0) * f[3];
                T x5 = (1.0 / 2.0) * f[6];
                e += dx *
                     ((1.0 / 2.0) * lambda * pow(x0 + x1 + x2 - 3.0 / 2.0, 2) +
                      mu * (pow(x0 - 1.0 / 2.0, 2) + pow(x1 - 1.0 / 2.0, 2) + pow(x2 - 1.0 / 2.0, 2) +
                            2 * pow((1.0 / 2.0) * f[1] * f[2] + (1.0 / 2.0) * f[4] * f[5] + (1.0 / 2.0) * f[7] * f[8],
                                    2) +
                            2 * pow(f[1] * x3 + f[4] * x4 + f[7] * x5, 2) +
                            2 * pow(f[2] * x3 + f[5] * x4 + f[8] * x5, 2)));
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
                T x1 = pow(f[3], 2);
                T x2 = pow(f[6], 2);
                T x3 = (1.0 / 2.0) * x0 + (1.0 / 2.0) * x1 + (1.0 / 2.0) * x2;
                T x4 = pow(f[1], 2);
                T x5 = pow(f[4], 2);
                T x6 = pow(f[7], 2);
                T x7 = (1.0 / 2.0) * x4 + (1.0 / 2.0) * x5 + (1.0 / 2.0) * x6;
                T x8 = pow(f[2], 2);
                T x9 = pow(f[5], 2);
                T x10 = pow(f[8], 2);
                T x11 = (1.0 / 2.0) * x10 + (1.0 / 2.0) * x8 + (1.0 / 2.0) * x9;
                T x12 = x11 + x3 + x7 - 3.0 / 2.0;
                T x13 = f[0] * f[1];
                T x14 = f[3] * f[4];
                T x15 = f[6] * f[7];
                T x16 = (1.0 / 2.0) * x13 + (1.0 / 2.0) * x14 + (1.0 / 2.0) * x15;
                T x17 = f[0] * f[2];
                T x18 = f[3] * f[5];
                T x19 = f[6] * f[8];
                T x20 = (1.0 / 2.0) * x17 + (1.0 / 2.0) * x18 + (1.0 / 2.0) * x19;
                T x21 = f[1] * f[2];
                T x22 = f[4] * f[5];
                T x23 = f[7] * f[8];
                T x24 = (1.0 / 2.0) * x21 + (1.0 / 2.0) * x22 + (1.0 / 2.0) * x23;
                T x25 = x3 - 1.0 / 2.0;
                T x26 = x7 - 1.0 / 2.0;
                T x27 = x11 - 1.0 / 2.0;
                T x28 = lambda * x12;
                T x29 = 2 * f[1];
                T x30 = 2 * f[2];
                T x31 = 2 * f[0];
                T x32 = 2 * f[4];
                T x33 = 2 * f[5];
                T x34 = 2 * f[3];
                T x35 = 2 * f[7];
                T x36 = 2 * f[8];
                T x37 = 2 * f[6];
                T x38 = lambda * x13 + mu * (2 * x13 + x14 + x15);
                T x39 = lambda * x17 + mu * (2 * x17 + x18 + x19);
                T x40 = x1 + x4;
                T x41 = x8 - 1;
                T x42 = x2 + x41;
                T x43 = lambda * x21 + mu * (2 * x21 + x22 + x23);
                T x44 = x0 + x5;
                T x45 = x9 - 1;
                T x46 = x0 + x10;
                T x47 = f[0] * mu;
                T x48 = f[3] * lambda;
                T x49 = f[1] * x48 + f[4] * x47;
                T x50 = f[2] * x48 + f[5] * x47;
                T x51 = f[0] * f[3];
                T x52 = f[1] * f[4];
                T x53 = f[2] * f[5];
                T x54 = grad_test[0] * (lambda * x51 + mu * (2 * x51 + x52 + x53));
                T x55 = f[4] * lambda;
                T x56 = f[1] * mu;
                T x57 = f[0] * x55 + f[3] * x56;
                T x58 = f[2] * x55 + f[5] * x56;
                T x59 = grad_test[1] * (lambda * x52 + mu * (x51 + 2 * x52 + x53));
                T x60 = f[5] * lambda;
                T x61 = f[2] * mu;
                T x62 = f[0] * x60 + f[3] * x61;
                T x63 = f[1] * x60 + f[4] * x61;
                T x64 = grad_test[2] * (lambda * x53 + mu * (x51 + x52 + 2 * x53));
                T x65 = f[6] * lambda;
                T x66 = f[1] * x65 + f[7] * x47;
                T x67 = f[2] * x65 + f[8] * x47;
                T x68 = f[0] * f[6];
                T x69 = f[1] * f[7];
                T x70 = f[2] * f[8];
                T x71 = grad_test[0] * (lambda * x68 + mu * (2 * x68 + x69 + x70));
                T x72 = f[7] * lambda;
                T x73 = f[0] * x72 + f[6] * x56;
                T x74 = f[2] * x72 + f[8] * x56;
                T x75 = grad_test[1] * (lambda * x69 + mu * (x68 + 2 * x69 + x70));
                T x76 = f[8] * lambda;
                T x77 = f[0] * x76 + f[6] * x61;
                T x78 = f[1] * x76 + f[7] * x61;
                T x79 = grad_test[2] * (lambda * x70 + mu * (x68 + x69 + 2 * x70));
                T x80 = lambda * x14 + mu * (x13 + 2 * x14 + x15);
                T x81 = lambda * x18 + mu * (x17 + 2 * x18 + x19);
                T x82 = lambda * x22 + mu * (x21 + 2 * x22 + x23);
                T x83 = x10 + x5;
                T x84 = f[3] * mu;
                T x85 = f[6] * x55 + f[7] * x84;
                T x86 = f[6] * x60 + f[8] * x84;
                T x87 = f[3] * f[6];
                T x88 = f[4] * f[7];
                T x89 = f[5] * f[8];
                T x90 = grad_test[0] * (lambda * x87 + mu * (2 * x87 + x88 + x89));
                T x91 = f[4] * mu;
                T x92 = f[6] * x91 + f[7] * x48;
                T x93 = f[7] * x60 + f[8] * x91;
                T x94 = grad_test[1] * (lambda * x88 + mu * (x87 + 2 * x88 + x89));
                T x95 = f[5] * mu;
                T x96 = f[6] * x95 + f[8] * x48;
                T x97 = f[7] * x95 + f[8] * x55;
                T x98 = grad_test[2] * (lambda * x89 + mu * (x87 + x88 + 2 * x89));
                T x99 = lambda * x15 + mu * (x13 + x14 + 2 * x15);
                T x100 = lambda * x19 + mu * (x17 + x18 + 2 * x19);
                T x101 = lambda * x23 + mu * (x21 + x22 + 2 * x23);
                e += dx *
                     ((1.0 / 2.0) * lambda * pow(x12, 2) + mu * (2 * pow(x16, 2) + 2 * pow(x20, 2) + 2 * pow(x24, 2) +
                                                                 pow(x25, 2) + pow(x26, 2) + pow(x27, 2)));
                lf[0] += dx * (grad_test[0] * (f[0] * x28 + mu * (x16 * x29 + x20 * x30 + x25 * x31)) +
                               grad_test[1] * (f[1] * x28 + mu * (x16 * x31 + x24 * x30 + x26 * x29)) +
                               grad_test[2] * (f[2] * x28 + mu * (x20 * x31 + x24 * x29 + x27 * x30)));
                lf[1] += dx * (grad_test[0] * (f[3] * x28 + mu * (x16 * x32 + x20 * x33 + x25 * x34)) +
                               grad_test[1] * (f[4] * x28 + mu * (x16 * x34 + x24 * x33 + x26 * x32)) +
                               grad_test[2] * (f[5] * x28 + mu * (x20 * x34 + x24 * x32 + x27 * x33)));
                lf[2] += dx * (grad_test[0] * (f[6] * x28 + mu * (x16 * x35 + x20 * x36 + x25 * x37)) +
                               grad_test[1] * (f[7] * x28 + mu * (x16 * x37 + x24 * x36 + x26 * x35)) +
                               grad_test[2] * (f[8] * x28 + mu * (x20 * x37 + x24 * x35 + x27 * x36)));
                bf[0] += dx * (grad_trial[0] * (grad_test[0] * (lambda * x0 + mu * (3 * x0 + x40 + x42) + x28) +
                                                grad_test[1] * x38 + grad_test[2] * x39) +
                               grad_trial[1] * (grad_test[0] * x38 +
                                                grad_test[1] * (lambda * x4 + mu * (3 * x4 + x41 + x44 + x6) + x28) +
                                                grad_test[2] * x43) +
                               grad_trial[2] * (grad_test[0] * x39 + grad_test[1] * x43 +
                                                grad_test[2] * (lambda * x8 + mu * (x4 + x45 + x46 + 3 * x8) + x28)));
                bf[1] += dx * (grad_trial[0] * (grad_test[1] * x49 + grad_test[2] * x50 + x54) +
                               grad_trial[1] * (grad_test[0] * x57 + grad_test[2] * x58 + x59) +
                               grad_trial[2] * (grad_test[0] * x62 + grad_test[1] * x63 + x64));
                bf[2] += dx * (grad_trial[0] * (grad_test[1] * x66 + grad_test[2] * x67 + x71) +
                               grad_trial[1] * (grad_test[0] * x73 + grad_test[2] * x74 + x75) +
                               grad_trial[2] * (grad_test[0] * x77 + grad_test[1] * x78 + x79));
                bf[3] += dx * (grad_trial[0] * (grad_test[1] * x57 + grad_test[2] * x62 + x54) +
                               grad_trial[1] * (grad_test[0] * x49 + grad_test[2] * x63 + x59) +
                               grad_trial[2] * (grad_test[0] * x50 + grad_test[1] * x58 + x64));
                bf[4] += dx * (grad_trial[0] * (grad_test[0] * (lambda * x1 + mu * (3 * x1 + x2 + x44 + x45) + x28) +
                                                grad_test[1] * x80 + grad_test[2] * x81) +
                               grad_trial[1] * (grad_test[0] * x80 +
                                                grad_test[1] * (lambda * x5 + mu * (x40 + x45 + 3 * x5 + x6) + x28) +
                                                grad_test[2] * x82) +
                               grad_trial[2] * (grad_test[0] * x81 + grad_test[1] * x82 +
                                                grad_test[2] * (lambda * x9 + mu * (x1 + x41 + x83 + 3 * x9) + x28)));
                bf[5] += dx * (grad_trial[0] * (grad_test[1] * x85 + grad_test[2] * x86 + x90) +
                               grad_trial[1] * (grad_test[0] * x92 + grad_test[2] * x93 + x94) +
                               grad_trial[2] * (grad_test[0] * x96 + grad_test[1] * x97 + x98));
                bf[6] += dx * (grad_trial[0] * (grad_test[1] * x73 + grad_test[2] * x77 + x71) +
                               grad_trial[1] * (grad_test[0] * x66 + grad_test[2] * x78 + x75) +
                               grad_trial[2] * (grad_test[0] * x67 + grad_test[1] * x74 + x79));
                bf[7] += dx * (grad_trial[0] * (grad_test[1] * x92 + grad_test[2] * x96 + x90) +
                               grad_trial[1] * (grad_test[0] * x85 + grad_test[2] * x97 + x94) +
                               grad_trial[2] * (grad_test[0] * x86 + grad_test[1] * x93 + x98));
                bf[8] += dx * (grad_trial[0] * (grad_test[0] * (lambda * x2 + mu * (x1 + 3 * x2 + x46 + x6 - 1) + x28) +
                                                grad_test[1] * x99 + grad_test[2] * x100) +
                               grad_trial[1] * (grad_test[0] * x99 +
                                                grad_test[1] * (lambda * x6 + mu * (x2 + x4 + 3 * x6 + x83 - 1) + x28) +
                                                grad_test[2] * x101) +
                               grad_trial[2] * (grad_test[0] * x100 + grad_test[1] * x101 +
                                                grad_test[2] * (lambda * x10 + mu * (3 * x10 + x42 + x6 + x9) + x28)));
            }

            T mu{1.0};
            T lambda{1.0};
        };
    }  // namespace kernels
}  // namespace utopia

#endif  // UTOPIA_TPL_HYPERELASTICITY_SaintVenantKirchoff_3_IMPL_hpp
