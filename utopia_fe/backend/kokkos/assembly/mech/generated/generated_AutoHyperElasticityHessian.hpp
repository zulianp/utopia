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
    T x17 = pow(x10 + x11 + x12 + x13 + x14 + x15 + x16 + x8 + x9 + 1, 2);
    T x18 = 11 * lmbda;
    T x19 = x17 * x18;
    T x20 = x19 * x7;
    T x21 = f[1] * x1 - x20 * x4;
    T x22 = f[3] * f[7];
    T x23 = f[4] * f[6];
    T x24 = x22 - x23;
    T x25 = f[2] * x1 + x20 * x24;
    T x26 = -x9;
    T x27 = -x10;
    T x28 = -x11;
    T x29 = -x12;
    T x30 = -x13;
    T x31 = -x14;
    T x32 = -x15;
    T x33 = -x16;
    T x34 = x17 + x27 + x28 + x29 + x30 + x31 + x32 + x33 - 1;
    T x35 = 8 * mu;
    T x36 = f[1] * x0;
    T x37 = x19 * x4;
    T x38 = f[2] * x36 - x24 * x37;
    T x39 = -x8;
    T x40 = x17 + x26 + x29 + x30 + x31 + x32 + x33 + x39 - 1;
    T x41 = (1.0 / 6.0) * dx / x17;
    T x42 = f[1] * f[8];
    T x43 = f[2] * f[7];
    T x44 = x42 - x43;
    T x45 = grad_test[0] * (f[3] * x1 - x20 * x44);
    T x46 = f[1] * f[3];
    T x47 = x18 * x44;
    T x48 = 6 * mu + x18 * (-f[0] * x5 + f[0] * x6 + f[1] * x2 - f[1] * x3 - f[2] * x22 + f[2] * x23 + 1);
    T x49 = f[8] * x48;
    T x50 = x0 * x46 + x17 * (x4 * x47 + x49);
    T x51 = f[2] * f[3];
    T x52 = f[7] * x48;
    T x53 = -x0 * x51 + x17 * (x24 * x47 + x52);
    T x54 = f[0] * f[8];
    T x55 = f[2] * f[6];
    T x56 = x54 - x55;
    T x57 = grad_test[1] * (f[4] * x36 - x37 * x56);
    T x58 = f[2] * f[4];
    T x59 = x18 * x56;
    T x60 = f[6] * x48;
    T x61 = x0 * x58 + x17 * (x24 * x59 + x60);
    T x62 = f[0] * f[4];
    T x63 = x0 * x62 + x17 * (-x49 + x59 * x7);
    T x64 = f[2] * x0;
    T x65 = f[0] * f[7];
    T x66 = f[1] * f[6];
    T x67 = x65 - x66;
    T x68 = x19 * x24;
    T x69 = grad_test[2] * (f[5] * x64 - x67 * x68);
    T x70 = f[0] * f[5];
    T x71 = x18 * x67;
    T x72 = x0 * x70 + x17 * (x52 - x7 * x71);
    T x73 = f[1] * f[5];
    T x74 = x0 * x73 + x17 * (x4 * x71 - x60);
    T x75 = -x58 + x73;
    T x76 = grad_test[0] * (f[6] * x1 + x20 * x75);
    T x77 = x18 * x75;
    T x78 = f[4] * x48;
    T x79 = x0 * x55 + x17 * (x24 * x77 + x78);
    T x80 = f[5] * x48;
    T x81 = -x0 * x66 + x17 * (x4 * x77 + x80);
    T x82 = -x51 + x70;
    T x83 = grad_test[1] * (f[7] * x36 + x37 * x82);
    T x84 = x18 * x82;
    T x85 = x0 * x65 + x17 * (-x7 * x84 + x80);
    T x86 = f[3] * x48;
    T x87 = -x0 * x43 + x17 * (x24 * x84 + x86);
    T x88 = -x46 + x62;
    T x89 = grad_test[2] * (f[8] * x64 + x68 * x88);
    T x90 = x18 * x88;
    T x91 = x0 * x42 + x17 * (-x4 * x90 + x86);
    T x92 = x0 * x54 + x17 * (x7 * x90 - x78);
    T x93 = f[3] * x0;
    T x94 = x19 * x44;
    T x95 = f[4] * x93 - x56 * x94;
    T x96 = f[5] * x93 + x67 * x94;
    T x97 = f[4] * x0;
    T x98 = x19 * x56;
    T x99 = f[5] * x97 - x67 * x98;
    T x100 = x17 + x26 + x27 + x28 + x31 + x32 + x33 + x39 - 1;
    T x101 = grad_test[0] * (f[6] * x93 - x75 * x94);
    T x102 = f[2] * x48;
    T x103 = x0 * x23 + x17 * (x102 + x59 * x75);
    T x104 = f[1] * x48;
    T x105 = -x0 * x3 + x17 * (x104 + x71 * x75);
    T x106 = grad_test[1] * (f[7] * x97 - x82 * x98);
    T x107 = f[0] * x48;
    T x108 = x0 * x6 + x17 * (x107 + x71 * x82);
    T x109 = x0 * x22 + x17 * (-x102 + x47 * x82);
    T x110 = f[8] * x0;
    T x111 = x19 * x88;
    T x112 = grad_test[2] * (f[5] * x110 - x111 * x67);
    T x113 = x0 * x2 + x17 * (x104 - x47 * x88);
    T x114 = x0 * x5 + x17 * (-x107 + x59 * x88);
    T x115 = f[6] * f[7] * x0 - x19 * x75 * x82;
    T x116 = f[6] * x110 + x111 * x75;
    T x117 = x17 + x26 + x27 + x28 + x29 + x30 + x33 + x39 - 1;
    T x118 = f[7] * x110 - x111 * x82;
    bf[0] += x41 * (grad_trial[0] * (grad_test[0] * (x19 * pow(x7, 2) + x35 * (x26 + x34 + x8)) + grad_test[1] * x21 +
                                     grad_test[2] * x25) +
                    grad_trial[1] * (grad_test[0] * x21 + grad_test[1] * (x19 * pow(x4, 2) + x35 * (x34 + x39 + x9)) +
                                     grad_test[2] * x38) +
                    grad_trial[2] * (grad_test[0] * x25 + grad_test[1] * x38 +
                                     grad_test[2] * (x19 * pow(x24, 2) + x35 * (x10 + x28 + x40))));
    bf[1] += x41 * (grad_trial[0] * (grad_test[1] * x50 - grad_test[2] * x53 + x45) +
                    grad_trial[1] * (grad_test[0] * x63 + grad_test[2] * x61 + x57) +
                    grad_trial[2] * (grad_test[0] * x72 + grad_test[1] * x74 + x69));
    bf[2] += x41 * (grad_trial[0] * (-grad_test[1] * x81 + grad_test[2] * x79 + x76) +
                    grad_trial[1] * (grad_test[0] * x85 - grad_test[2] * x87 + x83) +
                    grad_trial[2] * (grad_test[0] * x92 + grad_test[1] * x91 + x89));
    bf[3] += x41 * (grad_trial[0] * (grad_test[1] * x63 + grad_test[2] * x72 + x45) +
                    grad_trial[1] * (grad_test[0] * x50 + grad_test[2] * x74 + x57) +
                    grad_trial[2] * (-grad_test[0] * x53 + grad_test[1] * x61 + x69));
    bf[4] +=
        x41 *
        (grad_trial[0] *
             (grad_test[0] * (x19 * pow(x44, 2) + x35 * (x11 + x27 + x40)) + grad_test[1] * x95 + grad_test[2] * x96) +
         grad_trial[1] *
             (grad_test[0] * x95 + grad_test[1] * (x19 * pow(x56, 2) + x35 * (x100 + x12 + x30)) + grad_test[2] * x99) +
         grad_trial[2] *
             (grad_test[0] * x96 + grad_test[1] * x99 + grad_test[2] * (x19 * pow(x67, 2) + x35 * (x100 + x13 + x29))));
    bf[5] += x41 * (grad_trial[0] * (grad_test[1] * x103 - grad_test[2] * x105 + x101) +
                    grad_trial[1] * (grad_test[0] * x109 + grad_test[2] * x108 + x106) +
                    grad_trial[2] * (grad_test[0] * x113 + grad_test[1] * x114 + x112));
    bf[6] += x41 * (grad_trial[0] * (grad_test[1] * x85 + grad_test[2] * x92 + x76) +
                    grad_trial[1] * (-grad_test[0] * x81 + grad_test[2] * x91 + x83) +
                    grad_trial[2] * (grad_test[0] * x79 - grad_test[1] * x87 + x89));
    bf[7] += x41 * (grad_trial[0] * (grad_test[1] * x109 + grad_test[2] * x113 + x101) +
                    grad_trial[1] * (grad_test[0] * x103 + grad_test[2] * x114 + x106) +
                    grad_trial[2] * (-grad_test[0] * x105 + grad_test[1] * x108 + x112));
    bf[8] += x41 *
             (grad_trial[0] * (grad_test[0] * (x19 * pow(x75, 2) + x35 * (x117 + x14 + x32)) + grad_test[1] * x115 +
                               grad_test[2] * x116) +
              grad_trial[1] * (grad_test[0] * x115 + grad_test[1] * (x19 * pow(x82, 2) + x35 * (x117 + x15 + x31)) +
                               grad_test[2] * x118) +
              grad_trial[2] * (grad_test[0] * x116 + grad_test[1] * x118 +
                               grad_test[2] * (x19 * pow(x88, 2) +
                                               x35 * (x16 + x17 + x26 + x27 + x28 + x29 + x30 + x31 + x32 + x39 - 1))));
}

#undef UTOPIA_RESTRICT

#endif  // UTOPIA_TPL_ELASTICITY_HESSIAN_HPP
