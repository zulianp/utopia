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
    T x0 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) + pow(f[6], 2) +
           pow(f[7], 2) + pow(f[8], 2);
    T x1 = 8 * mu * x0;
    T x2 = f[3] * f[8];
    T x3 = f[5] * f[6];
    T x4 = f[5] * f[7];
    T x5 = f[4] * f[6];
    T x6 = f[4] * f[8];
    T x7 = f[3] * f[7];
    T x8 = x0 + 1;
    T x9 = x8 * (11 * lmbda * (f[0] * x4 - f[0] * x6 + f[1] * x2 - f[1] * x3 + f[2] * x5 - f[2] * x7 + 1) + 6 * mu);
    T x10 = (1.0 / 6.0) * dx / x8;
    lf[0] += x10 * (grad_test[0] * (f[0] * x1 - x9 * (-x4 + x6)) + grad_test[1] * (f[1] * x1 + x9 * (x2 - x3)) +
                    grad_test[2] * (f[2] * x1 - x9 * (-x5 + x7)));
    lf[1] += x10 * (grad_test[0] * (f[3] * x1 + x9 * (f[1] * f[8] - f[2] * f[7])) +
                    grad_test[1] * (f[4] * x1 - x9 * (f[0] * f[8] - f[2] * f[6])) +
                    grad_test[2] * (f[5] * x1 + x9 * (f[0] * f[7] - f[1] * f[6])));
    lf[2] += x10 * (grad_test[0] * (f[6] * x1 - x9 * (f[1] * f[5] - f[2] * f[4])) +
                    grad_test[1] * (f[7] * x1 + x9 * (f[0] * f[5] - f[2] * f[3])) +
                    grad_test[2] * (f[8] * x1 - x9 * (f[0] * f[4] - f[1] * f[3])));
}

#undef UTOPIA_RESTRICT

#endif  // UTOPIA_TPL_ELASTICITY_GRADIENT_HPP
