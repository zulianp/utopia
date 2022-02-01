#ifndef UTOPIA_TPL_ELASTICITY_HESSIAN_HPP
#define UTOPIA_TPL_ELASTICITY_HESSIAN_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

#ifndef AUTO_HYPER_ELASTICITY_DIM
#define AUTO_HYPER_ELASTICITY_DIM 3
#endif  // AUTO_HYPER_ELASTICITY_DIM

template <typename T>
UTOPIA_FUNCTION void elastic_material_hessian(const T mu,
                                              const T lmbda,
                                              const T *UTOPIA_RESTRICT f,
                                              const T *grad_test,
                                              const T *grad_trial,
                                              const T dx,
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
    T x17 = x10 + x11 + x12 + x13 + x14 + x15 + x16 + x8 + x9 + 1;
    T x18 = pow(x17, -2);
    T x19 = lmbda * mu;
    T x20 = 352 * x18 * x19;
    T x21 = f[0] * x20;
    T x22 = f[1] * x21 + x2 * x7;
    T x23 = (1.0 / 132.0) / lmbda;
    T x24 = grad_test[1] * x23;
    T x25 = f[4] * f[6];
    T x26 = f[3] * f[7];
    T x27 = x25 - x26;
    T x28 = f[2] * x21 + x27 * x7;
    T x29 = grad_test[2] * x23;
    T x30 = 4 * x18;
    T x31 = 2 - 2 / x17;
    T x32 = 88 * x19;
    T x33 = grad_test[0] * x23;
    T x34 = x2 * x6;
    T x35 = f[1] * x20;
    T x36 = f[2] * x35 + x27 * x34;
    T x37 = f[1] * f[8];
    T x38 = f[2] * f[7];
    T x39 = x37 - x38;
    T x40 = x33 * (f[3] * x21 + x39 * x7);
    T x41 = f[1] * f[3];
    T x42 = 22 * lmbda *
            (11 * lmbda * (f[0] * x3 - f[0] * x4 + f[1] * x0 - f[1] * x1 + f[2] * x25 - f[2] * x26 + 1) + 6 * mu);
    T x43 = f[8] * x42;
    T x44 = x20 * x41 + x34 * x39 + x43;
    T x45 = x27 * x6;
    T x46 = f[2] * f[3];
    T x47 = f[7] * x42;
    T x48 = x20 * x46 + x39 * x45 - x47;
    T x49 = f[2] * f[6];
    T x50 = f[0] * f[8];
    T x51 = x49 - x50;
    T x52 = x24 * (f[4] * x35 + x34 * x51);
    T x53 = f[0] * f[4];
    T x54 = x20 * x53 - x43 + x51 * x7;
    T x55 = f[2] * f[4];
    T x56 = f[6] * x42;
    T x57 = x20 * x55 + x45 * x51 + x56;
    T x58 = f[0] * f[7];
    T x59 = f[1] * f[6];
    T x60 = x58 - x59;
    T x61 = f[2] * x20;
    T x62 = x29 * (f[5] * x61 + x45 * x60);
    T x63 = f[0] * f[5];
    T x64 = x20 * x63 + x47 + x60 * x7;
    T x65 = f[1] * f[5];
    T x66 = x20 * x65 + x34 * x60 - x56;
    T x67 = x55 - x65;
    T x68 = x33 * (f[6] * x21 + x67 * x7);
    T x69 = f[5] * x42;
    T x70 = x20 * x59 + x34 * x67 - x69;
    T x71 = f[4] * x42;
    T x72 = x20 * x49 + x45 * x67 + x71;
    T x73 = -x46 + x63;
    T x74 = x24 * (f[7] * x35 + x34 * x73);
    T x75 = x20 * x58 + x69 + x7 * x73;
    T x76 = f[3] * x42;
    T x77 = x20 * x38 + x45 * x73 - x76;
    T x78 = x41 - x53;
    T x79 = x29 * (f[8] * x61 + x45 * x78);
    T x80 = x20 * x50 + x7 * x78 - x71;
    T x81 = x20 * x37 + x34 * x78 + x76;
    T x82 = x39 * x6;
    T x83 = f[3] * x20;
    T x84 = f[4] * x83 + x51 * x82;
    T x85 = f[5] * x83 + x60 * x82;
    T x86 = x51 * x6;
    T x87 = f[4] * x20;
    T x88 = f[5] * x87 + x60 * x86;
    T x89 = x33 * (f[6] * x83 + x67 * x82);
    T x90 = f[2] * x42;
    T x91 = x20 * x25 + x67 * x86 + x90;
    T x92 = x6 * x60;
    T x93 = f[1] * x42;
    T x94 = x1 * x20 + x67 * x92 - x93;
    T x95 = x24 * (f[7] * x87 + x73 * x86);
    T x96 = x20 * x26 + x73 * x82 - x90;
    T x97 = f[0] * x42;
    T x98 = x20 * x3 + x73 * x92 + x97;
    T x99 = f[8] * x20;
    T x100 = x29 * (f[5] * x99 + x78 * x92);
    T x101 = x0 * x20 + x78 * x82 + x93;
    T x102 = x20 * x4 + x78 * x86 - x97;
    T x103 = x6 * x67;
    T x104 = f[6] * f[7] * x20 + x103 * x73;
    T x105 = f[6] * x99 + x103 * x78;
    T x106 = f[7] * x99 + x6 * x73 * x78;
    bf[0] += dx * (grad_trial[0] * (x22 * x24 + x28 * x29 + x33 * (x32 * (x30 * x8 + x31) + pow(x5, 2) * x6)) +
                   grad_trial[1] * (x22 * x33 + x24 * (pow(x2, 2) * x6 + x32 * (x30 * x9 + x31)) + x29 * x36) +
                   grad_trial[2] * (x24 * x36 + x28 * x33 + x29 * (pow(x27, 2) * x6 + x32 * (x10 * x30 + x31))));
    bf[1] += dx * (grad_trial[0] * (x24 * x44 + x29 * x48 + x40) + grad_trial[1] * (x29 * x57 + x33 * x54 + x52) +
                   grad_trial[2] * (x24 * x66 + x33 * x64 + x62));
    bf[2] += dx * (grad_trial[0] * (x24 * x70 + x29 * x72 + x68) + grad_trial[1] * (x29 * x77 + x33 * x75 + x74) +
                   grad_trial[2] * (x24 * x81 + x33 * x80 + x79));
    bf[3] += dx * (grad_trial[0] * (x24 * x54 + x29 * x64 + x40) + grad_trial[1] * (x29 * x66 + x33 * x44 + x52) +
                   grad_trial[2] * (x24 * x57 + x33 * x48 + x62));
    bf[4] += dx * (grad_trial[0] * (x24 * x84 + x29 * x85 + x33 * (x32 * (x11 * x30 + x31) + pow(x39, 2) * x6)) +
                   grad_trial[1] * (x24 * (x32 * (x12 * x30 + x31) + pow(x51, 2) * x6) + x29 * x88 + x33 * x84) +
                   grad_trial[2] * (x24 * x88 + x29 * (x32 * (x13 * x30 + x31) + x6 * pow(x60, 2)) + x33 * x85));
    bf[5] += dx * (grad_trial[0] * (x24 * x91 + x29 * x94 + x89) + grad_trial[1] * (x29 * x98 + x33 * x96 + x95) +
                   grad_trial[2] * (x100 + x101 * x33 + x102 * x24));
    bf[6] += dx * (grad_trial[0] * (x24 * x75 + x29 * x80 + x68) + grad_trial[1] * (x29 * x81 + x33 * x70 + x74) +
                   grad_trial[2] * (x24 * x77 + x33 * x72 + x79));
    bf[7] += dx * (grad_trial[0] * (x101 * x29 + x24 * x96 + x89) + grad_trial[1] * (x102 * x29 + x33 * x91 + x95) +
                   grad_trial[2] * (x100 + x24 * x98 + x33 * x94));
    bf[8] += dx * (grad_trial[0] * (x104 * x24 + x105 * x29 + x33 * (x32 * (x14 * x30 + x31) + x6 * pow(x67, 2))) +
                   grad_trial[1] * (x104 * x33 + x106 * x29 + x24 * (x32 * (x15 * x30 + x31) + x6 * pow(x73, 2))) +
                   grad_trial[2] * (x105 * x33 + x106 * x24 + x29 * (x32 * (x16 * x30 + x31) + x6 * pow(x78, 2))));
}

#undef UTOPIA_RESTRICT

#endif  // UTOPIA_TPL_ELASTICITY_HESSIAN_HPP
