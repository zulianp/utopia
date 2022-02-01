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
    T x0 = f[3] * f[8];
    T x1 = f[5] * f[6];
    T x2 = x0 - x1;
    T x3 = f[5] * f[7];
    T x4 = f[4] * f[8];
    T x5 = x3 - x4;
    T x6 = 242 * pow(lmbda, 2);
    T x7 = x5 * x6;
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
    T x19 = pow(x18, -2);
    T x20 = lmbda * mu;
    T x21 = 352 * x19 * x20;
    T x22 = f[0] * x21;
    T x23 = f[1] * x22 + x2 * x7;
    T x24 = (1.0 / 132.0) / lmbda;
    T x25 = grad_test[1] * x24;
    T x26 = f[4] * f[6];
    T x27 = f[3] * f[7];
    T x28 = x26 - x27;
    T x29 = f[2] * x22 + x28 * x7;
    T x30 = grad_test[2] * x24;
    T x31 = 4 * x19;
    T x32 = 2 / x18;
    T x33 = 2 - x32;
    T x34 = 88 * x20;
    T x35 = grad_test[0] * x24;
    T x36 = x2 * x6;
    T x37 = f[1] * x21;
    T x38 = f[2] * x37 + x28 * x36;
    T x39 = f[1] * f[8];
    T x40 = f[2] * f[7];
    T x41 = x39 - x40;
    T x42 = x35 * (f[3] * x22 + x41 * x7);
    T x43 = f[1] * f[3];
    T x44 = 11 * lmbda * (f[0] * x3 - f[0] * x4 + f[1] * x0 - f[1] * x1 + f[2] * x26 - f[2] * x27 + 1) + 6 * mu;
    T x45 = 22 * lmbda * x44;
    T x46 = f[8] * x45;
    T x47 = x21 * x43 + x36 * x41 + x46;
    T x48 = x28 * x6;
    T x49 = f[2] * f[3];
    T x50 = f[7] * x45;
    T x51 = x21 * x49 + x41 * x48 - x50;
    T x52 = f[2] * f[6];
    T x53 = f[0] * f[8];
    T x54 = x52 - x53;
    T x55 = x25 * (f[4] * x37 + x36 * x54);
    T x56 = f[0] * f[4];
    T x57 = x21 * x56 - x46 + x54 * x7;
    T x58 = f[2] * f[4];
    T x59 = f[6] * x45;
    T x60 = x21 * x58 + x48 * x54 + x59;
    T x61 = f[0] * f[7];
    T x62 = f[1] * f[6];
    T x63 = x61 - x62;
    T x64 = f[2] * x21;
    T x65 = x30 * (f[5] * x64 + x48 * x63);
    T x66 = f[0] * f[5];
    T x67 = x21 * x66 + x50 + x63 * x7;
    T x68 = f[1] * f[5];
    T x69 = x21 * x68 + x36 * x63 - x59;
    T x70 = x58 - x68;
    T x71 = x35 * (f[6] * x22 + x7 * x70);
    T x72 = f[5] * x45;
    T x73 = x21 * x62 + x36 * x70 - x72;
    T x74 = f[4] * x45;
    T x75 = x21 * x52 + x48 * x70 + x74;
    T x76 = -x49 + x66;
    T x77 = x25 * (f[7] * x37 + x36 * x76);
    T x78 = x21 * x61 + x7 * x76 + x72;
    T x79 = f[3] * x45;
    T x80 = x21 * x40 + x48 * x76 - x79;
    T x81 = x43 - x56;
    T x82 = x30 * (f[8] * x64 + x48 * x81);
    T x83 = x21 * x53 + x7 * x81 - x74;
    T x84 = x21 * x39 + x36 * x81 + x79;
    T x85 = x41 * x6;
    T x86 = f[3] * x21;
    T x87 = f[4] * x86 + x54 * x85;
    T x88 = f[5] * x86 + x63 * x85;
    T x89 = x54 * x6;
    T x90 = f[4] * x21;
    T x91 = f[5] * x90 + x63 * x89;
    T x92 = x35 * (f[6] * x86 + x70 * x85);
    T x93 = f[2] * x45;
    T x94 = x21 * x26 + x70 * x89 + x93;
    T x95 = x6 * x63;
    T x96 = f[1] * x45;
    T x97 = x1 * x21 + x70 * x95 - x96;
    T x98 = x25 * (f[7] * x90 + x76 * x89);
    T x99 = x21 * x27 + x76 * x85 - x93;
    T x100 = f[0] * x45;
    T x101 = x100 + x21 * x3 + x76 * x95;
    T x102 = f[8] * x21;
    T x103 = x30 * (f[5] * x102 + x81 * x95);
    T x104 = x0 * x21 + x81 * x85 + x96;
    T x105 = -x100 + x21 * x4 + x81 * x89;
    T x106 = x6 * x70;
    T x107 = f[6] * f[7] * x21 + x106 * x76;
    T x108 = f[6] * x102 + x106 * x81;
    T x109 = f[7] * x102 + x6 * x76 * x81;
    bf[0] += dx * (grad_trial[0] * (x23 * x25 + x29 * x30 + x35 * (x34 * (x31 * x8 + x33) + pow(x5, 2) * x6)) +
                   grad_trial[1] * (x23 * x35 + x25 * (pow(x2, 2) * x6 + x34 * (x31 * x9 + x33)) + x30 * x38) +
                   grad_trial[2] * (x25 * x38 + x29 * x35 + x30 * (pow(x28, 2) * x6 + x34 * (x10 * x31 + x33))));
    bf[1] += dx * (grad_trial[0] * (x25 * x47 + x30 * x51 + x42) + grad_trial[1] * (x30 * x60 + x35 * x57 + x55) +
                   grad_trial[2] * (x25 * x69 + x35 * x67 + x65));
    bf[2] += dx * (grad_trial[0] * (x25 * x73 + x30 * x75 + x71) + grad_trial[1] * (x30 * x80 + x35 * x78 + x77) +
                   grad_trial[2] * (x25 * x84 + x35 * x83 + x82));
    bf[3] += dx * (grad_trial[0] * (x25 * x57 + x30 * x67 + x42) + grad_trial[1] * (x30 * x69 + x35 * x47 + x55) +
                   grad_trial[2] * (x25 * x60 + x35 * x51 + x65));
    bf[4] += dx * (grad_trial[0] * (x25 * x87 + x30 * x88 + x35 * (x34 * (x11 * x31 + x33) + pow(x41, 2) * x6)) +
                   grad_trial[1] * (x25 * (x34 * (x12 * x31 + x33) + pow(x54, 2) * x6) + x30 * x91 + x35 * x87) +
                   grad_trial[2] * (x25 * x91 + x30 * (x34 * (x13 * x31 + x33) + x6 * pow(x63, 2)) + x35 * x88));
    bf[5] += dx * (grad_trial[0] * (x25 * x94 + x30 * x97 + x92) + grad_trial[1] * (x101 * x30 + x35 * x99 + x98) +
                   grad_trial[2] * (x103 + x104 * x35 + x105 * x25));
    bf[6] += dx * (grad_trial[0] * (x25 * x78 + x30 * x83 + x71) + grad_trial[1] * (x30 * x84 + x35 * x73 + x77) +
                   grad_trial[2] * (x25 * x80 + x35 * x75 + x82));
    bf[7] += dx * (grad_trial[0] * (x104 * x30 + x25 * x99 + x92) + grad_trial[1] * (x105 * x30 + x35 * x94 + x98) +
                   grad_trial[2] * (x101 * x25 + x103 + x35 * x97));
    bf[8] += dx * (grad_trial[0] * (x107 * x25 + x108 * x30 + x35 * (x34 * (x14 * x31 + x33) + x6 * pow(x70, 2))) +
                   grad_trial[1] * (x107 * x35 + x109 * x30 + x25 * (x34 * (x15 * x31 + x33) + x6 * pow(x76, 2))) +
                   grad_trial[2] * (x108 * x35 + x109 * x25 + x30 * (x34 * (x16 * x31 + x33) + x6 * pow(x81, 2))));
    lf[0] +=
        dx * (x25 * (x2 * x45 + x34 * (-f[1] * x32 + 2 * f[1])) + x30 * (x28 * x45 + x34 * (-f[2] * x32 + 2 * f[2])) +
              x35 * (x34 * (-f[0] * x32 + 2 * f[0]) + x45 * x5));
    lf[1] +=
        dx * (x25 * (x34 * (-f[4] * x32 + 2 * f[4]) + x45 * x54) + x30 * (x34 * (-f[5] * x32 + 2 * f[5]) + x45 * x63) +
              x35 * (x34 * (-f[3] * x32 + 2 * f[3]) + x41 * x45));
    lf[2] +=
        dx * (x25 * (x34 * (-f[7] * x32 + 2 * f[7]) + x45 * x76) + x30 * (x34 * (-f[8] * x32 + 2 * f[8]) + x45 * x81) +
              x35 * (x34 * (-f[6] * x32 + 2 * f[6]) + x45 * x70));
    e += dx * x24 * (x34 * (x17 - log(x18) - 3) + pow(x44, 2));
}

#undef UTOPIA_RESTRICT

#endif  // UTOPIA_TPL_ELASTICITY_HPP
