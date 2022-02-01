#ifndef UTOPIA_TPL_ELASTICITY_HPP
#define UTOPIA_TPL_ELASTICITY_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

#ifndef AUTO_HYPER_ELASTICITY_DIM
#define AUTO_HYPER_ELASTICITY_DIM 3
#endif  // AUTO_HYPER_ELASTICITY_DIM

template <typename T>
UTOPIA_FUNCTION void elastic_material(const T mu,
                                      const T lmbda,
                                      const T *UTOPIA_RESTRICT f,
                                      const T *grad_test,
                                      const T *grad_trial,
                                      const T dx,
                                      T &e,
                                      T *UTOPIA_RESTRICT lf,
                                      T *UTOPIA_RESTRICT bf) {
    using namespace utopia::device;
    // Automatically generated
    T x0 = 16 * mu;
    T x1 = f[0] * x0;
    T x2 = f[3] * f[8];
    T x3 = f[5] * f[6];
    T x4 = x2 - x3;
    T x5 = f[4] * f[8];
    T x6 = f[5] * f[7];
    T x7 = x5 - x6;
    T x8 = pow(f[0], 2);
    T x9 = pow(f[1], 2);
    T x10 = pow(f[2], 2);
    T x11 = pow(f[3], 2);
    T x12 = pow(f[4], 2);
    T x13 = pow(f[5], 2);
    T x14 = pow(f[6], 2);
    T x15 = pow(f[7], 2);
    T x16 = pow(f[8], 2);
    T x17 = x10 + x11 + x12 + x13 + x14 + x15 + x16 + x8 + x9;
    T x18 = x17 + 1;
    T x19 = pow(x18, 2);
    T x20 = 11 * lmbda;
    T x21 = x19 * x20;
    T x22 = x21 * x7;
    T x23 = f[1] * x1 - x22 * x4;
    T x24 = f[3] * f[7];
    T x25 = f[4] * f[6];
    T x26 = x24 - x25;
    T x27 = f[2] * x1 + x22 * x26;
    T x28 = -x9;
    T x29 = -x10;
    T x30 = -x11;
    T x31 = -x12;
    T x32 = -x13;
    T x33 = -x14;
    T x34 = -x15;
    T x35 = -x16;
    T x36 = x19 + x29 + x30 + x31 + x32 + x33 + x34 + x35 - 1;
    T x37 = 8 * mu;
    T x38 = f[1] * x0;
    T x39 = x21 * x4;
    T x40 = f[2] * x38 - x26 * x39;
    T x41 = -x8;
    T x42 = x19 + x28 + x31 + x32 + x33 + x34 + x35 + x41 - 1;
    T x43 = (1.0 / 6.0) * dx;
    T x44 = x43 / x19;
    T x45 = f[1] * f[8];
    T x46 = f[2] * f[7];
    T x47 = x45 - x46;
    T x48 = grad_test[0] * (f[3] * x1 - x22 * x47);
    T x49 = f[1] * f[3];
    T x50 = x20 * x47;
    T x51 = 6 * mu + x20 * (-f[0] * x5 + f[0] * x6 + f[1] * x2 - f[1] * x3 - f[2] * x24 + f[2] * x25 + 1);
    T x52 = f[8] * x51;
    T x53 = x0 * x49 + x19 * (x4 * x50 + x52);
    T x54 = f[2] * f[3];
    T x55 = f[7] * x51;
    T x56 = -x0 * x54 + x19 * (x26 * x50 + x55);
    T x57 = f[0] * f[8];
    T x58 = f[2] * f[6];
    T x59 = x57 - x58;
    T x60 = grad_test[1] * (f[4] * x38 - x39 * x59);
    T x61 = f[2] * f[4];
    T x62 = x20 * x59;
    T x63 = f[6] * x51;
    T x64 = x0 * x61 + x19 * (x26 * x62 + x63);
    T x65 = f[0] * f[4];
    T x66 = x0 * x65 + x19 * (-x52 + x62 * x7);
    T x67 = f[2] * x0;
    T x68 = f[0] * f[7];
    T x69 = f[1] * f[6];
    T x70 = x68 - x69;
    T x71 = x21 * x26;
    T x72 = grad_test[2] * (f[5] * x67 - x70 * x71);
    T x73 = f[0] * f[5];
    T x74 = x20 * x70;
    T x75 = x0 * x73 + x19 * (x55 - x7 * x74);
    T x76 = f[1] * f[5];
    T x77 = x0 * x76 + x19 * (x4 * x74 - x63);
    T x78 = -x61 + x76;
    T x79 = grad_test[0] * (f[6] * x1 + x22 * x78);
    T x80 = x20 * x78;
    T x81 = f[4] * x51;
    T x82 = x0 * x58 + x19 * (x26 * x80 + x81);
    T x83 = f[5] * x51;
    T x84 = -x0 * x69 + x19 * (x4 * x80 + x83);
    T x85 = -x54 + x73;
    T x86 = grad_test[1] * (f[7] * x38 + x39 * x85);
    T x87 = x20 * x85;
    T x88 = x0 * x68 + x19 * (-x7 * x87 + x83);
    T x89 = f[3] * x51;
    T x90 = -x0 * x46 + x19 * (x26 * x87 + x89);
    T x91 = -x49 + x65;
    T x92 = grad_test[2] * (f[8] * x67 + x71 * x91);
    T x93 = x20 * x91;
    T x94 = x0 * x45 + x19 * (-x4 * x93 + x89);
    T x95 = x0 * x57 + x19 * (x7 * x93 - x81);
    T x96 = f[3] * x0;
    T x97 = x21 * x47;
    T x98 = f[4] * x96 - x59 * x97;
    T x99 = f[5] * x96 + x70 * x97;
    T x100 = f[4] * x0;
    T x101 = x21 * x59;
    T x102 = f[5] * x100 - x101 * x70;
    T x103 = x19 + x28 + x29 + x30 + x33 + x34 + x35 + x41 - 1;
    T x104 = grad_test[0] * (f[6] * x96 - x78 * x97);
    T x105 = f[2] * x51;
    T x106 = x0 * x25 + x19 * (x105 + x62 * x78);
    T x107 = f[1] * x51;
    T x108 = -x0 * x3 + x19 * (x107 + x74 * x78);
    T x109 = grad_test[1] * (f[7] * x100 - x101 * x85);
    T x110 = f[0] * x51;
    T x111 = x0 * x6 + x19 * (x110 + x74 * x85);
    T x112 = x0 * x24 + x19 * (-x105 + x50 * x85);
    T x113 = f[8] * x0;
    T x114 = x21 * x91;
    T x115 = grad_test[2] * (f[5] * x113 - x114 * x70);
    T x116 = x0 * x2 + x19 * (x107 - x50 * x91);
    T x117 = x0 * x5 + x19 * (-x110 + x62 * x91);
    T x118 = f[6] * f[7] * x0 - x21 * x78 * x85;
    T x119 = f[6] * x113 + x114 * x78;
    T x120 = x19 + x28 + x29 + x30 + x31 + x32 + x35 + x41 - 1;
    T x121 = f[7] * x113 - x114 * x85;
    T x122 = x17 * x37;
    T x123 = x18 * x51;
    T x124 = x43 / x18;
    bf[0] += x44 * (grad_trial[0] * (grad_test[0] * (x21 * pow(x7, 2) + x37 * (x28 + x36 + x8)) + grad_test[1] * x23 +
                                     grad_test[2] * x27) +
                    grad_trial[1] * (grad_test[0] * x23 + grad_test[1] * (x21 * pow(x4, 2) + x37 * (x36 + x41 + x9)) +
                                     grad_test[2] * x40) +
                    grad_trial[2] * (grad_test[0] * x27 + grad_test[1] * x40 +
                                     grad_test[2] * (x21 * pow(x26, 2) + x37 * (x10 + x30 + x42))));
    bf[1] += x44 * (grad_trial[0] * (grad_test[1] * x53 - grad_test[2] * x56 + x48) +
                    grad_trial[1] * (grad_test[0] * x66 + grad_test[2] * x64 + x60) +
                    grad_trial[2] * (grad_test[0] * x75 + grad_test[1] * x77 + x72));
    bf[2] += x44 * (grad_trial[0] * (-grad_test[1] * x84 + grad_test[2] * x82 + x79) +
                    grad_trial[1] * (grad_test[0] * x88 - grad_test[2] * x90 + x86) +
                    grad_trial[2] * (grad_test[0] * x95 + grad_test[1] * x94 + x92));
    bf[3] += x44 * (grad_trial[0] * (grad_test[1] * x66 + grad_test[2] * x75 + x48) +
                    grad_trial[1] * (grad_test[0] * x53 + grad_test[2] * x77 + x60) +
                    grad_trial[2] * (-grad_test[0] * x56 + grad_test[1] * x64 + x72));
    bf[4] +=
        x44 * (grad_trial[0] * (grad_test[0] * (x21 * pow(x47, 2) + x37 * (x11 + x29 + x42)) + grad_test[1] * x98 +
                                grad_test[2] * x99) +
               grad_trial[1] * (grad_test[0] * x98 + grad_test[1] * (x21 * pow(x59, 2) + x37 * (x103 + x12 + x32)) +
                                grad_test[2] * x102) +
               grad_trial[2] * (grad_test[0] * x99 + grad_test[1] * x102 +
                                grad_test[2] * (x21 * pow(x70, 2) + x37 * (x103 + x13 + x31))));
    bf[5] += x44 * (grad_trial[0] * (grad_test[1] * x106 - grad_test[2] * x108 + x104) +
                    grad_trial[1] * (grad_test[0] * x112 + grad_test[2] * x111 + x109) +
                    grad_trial[2] * (grad_test[0] * x116 + grad_test[1] * x117 + x115));
    bf[6] += x44 * (grad_trial[0] * (grad_test[1] * x88 + grad_test[2] * x95 + x79) +
                    grad_trial[1] * (-grad_test[0] * x84 + grad_test[2] * x94 + x86) +
                    grad_trial[2] * (grad_test[0] * x82 - grad_test[1] * x90 + x92));
    bf[7] += x44 * (grad_trial[0] * (grad_test[1] * x112 + grad_test[2] * x116 + x104) +
                    grad_trial[1] * (grad_test[0] * x106 + grad_test[2] * x117 + x109) +
                    grad_trial[2] * (-grad_test[0] * x108 + grad_test[1] * x111 + x115));
    bf[8] += x44 *
             (grad_trial[0] * (grad_test[0] * (x21 * pow(x78, 2) + x37 * (x120 + x14 + x34)) + grad_test[1] * x118 +
                               grad_test[2] * x119) +
              grad_trial[1] * (grad_test[0] * x118 + grad_test[1] * (x21 * pow(x85, 2) + x37 * (x120 + x15 + x33)) +
                               grad_test[2] * x121) +
              grad_trial[2] * (grad_test[0] * x119 + grad_test[1] * x121 +
                               grad_test[2] * (x21 * pow(x91, 2) +
                                               x37 * (x16 + x19 + x28 + x29 + x30 + x31 + x32 + x33 + x34 + x41 - 1))));
    lf[0] += x124 * (grad_test[0] * (f[0] * x122 - x123 * x7) + grad_test[1] * (f[1] * x122 + x123 * x4) +
                     grad_test[2] * (f[2] * x122 - x123 * x26));
    lf[1] += x124 * (grad_test[0] * (f[3] * x122 + x123 * x47) + grad_test[1] * (f[4] * x122 - x123 * x59) +
                     grad_test[2] * (f[5] * x122 + x123 * x70));
    lf[2] += x124 * (grad_test[0] * (f[6] * x122 - x123 * x78) + grad_test[1] * (f[7] * x122 + x123 * x85) +
                     grad_test[2] * (f[8] * x122 - x123 * x91));
    e += (1.0 / 132.0) * dx * (88 * lmbda * mu * (x17 - log(x18) - 3) + pow(x51, 2)) / lmbda;
}

#undef UTOPIA_RESTRICT

#endif  // UTOPIA_TPL_ELASTICITY_HPP
