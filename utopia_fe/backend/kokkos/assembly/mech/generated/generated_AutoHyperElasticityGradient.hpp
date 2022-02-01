#ifndef UTOPIA_TPL_ELASTICITY_GRADIENT_HPP
#define UTOPIA_TPL_ELASTICITY_GRADIENT_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

#ifndef AUTO_HYPER_ELASTICITY_DIM
#define AUTO_HYPER_ELASTICITY_DIM 3
#endif  // AUTO_HYPER_ELASTICITY_DIM

template <typename T>
UTOPIA_FUNCTION void elastic_material_gradient(const T mu,
                                               const T lmbda,
                                               const T *UTOPIA_RESTRICT f,
                                               const T *grad_test,
                                               const T dx,
                                               T *UTOPIA_RESTRICT lf) {
    using namespace utopia::device;
    // Automatically generated
    T x0 = 2 * f[0];
    T x1 = 1.0 / (pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) +
                  pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2) + 1);
    T x2 = 88 * lmbda * mu;
    T x3 = f[5] * f[7];
    T x4 = f[4] * f[8];
    T x5 = f[3] * f[8];
    T x6 = f[4] * f[6];
    T x7 = f[5] * f[6];
    T x8 = f[3] * f[7];
    T x9 = 22 * lmbda *
           (11 * lmbda * (f[0] * x3 - f[0] * x4 + f[1] * x5 - f[1] * x7 + f[2] * x6 - f[2] * x8 + 1) + 6 * mu);
    T x10 = (1.0 / 132.0) / lmbda;
    T x11 = grad_test[0] * x10;
    T x12 = 2 * f[1];
    T x13 = grad_test[1] * x10;
    T x14 = 2 * f[2];
    T x15 = grad_test[2] * x10;
    T x16 = 2 * f[3];
    T x17 = 2 * f[4];
    T x18 = 2 * f[5];
    T x19 = 2 * f[6];
    T x20 = 2 * f[7];
    T x21 = 2 * f[8];
    lf[0] += dx * (x11 * (x2 * (-x0 * x1 + x0) + x9 * (x3 - x4)) + x13 * (x2 * (-x1 * x12 + x12) + x9 * (x5 - x7)) +
                   x15 * (x2 * (-x1 * x14 + x14) + x9 * (x6 - x8)));
    lf[1] += dx * (x11 * (x2 * (-x1 * x16 + x16) + x9 * (f[1] * f[8] - f[2] * f[7])) +
                   x13 * (x2 * (-x1 * x17 + x17) + x9 * (-f[0] * f[8] + f[2] * f[6])) +
                   x15 * (x2 * (-x1 * x18 + x18) + x9 * (f[0] * f[7] - f[1] * f[6])));
    lf[2] += dx * (x11 * (x2 * (-x1 * x19 + x19) + x9 * (-f[1] * f[5] + f[2] * f[4])) +
                   x13 * (x2 * (-x1 * x20 + x20) + x9 * (f[0] * f[5] - f[2] * f[3])) +
                   x15 * (x2 * (-x1 * x21 + x21) + x9 * (-f[0] * f[4] + f[1] * f[3])));
}

#undef UTOPIA_RESTRICT

#endif  // UTOPIA_TPL_ELASTICITY_GRADIENT_HPP
