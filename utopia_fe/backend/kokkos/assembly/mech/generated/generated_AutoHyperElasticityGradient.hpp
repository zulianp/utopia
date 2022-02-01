#ifndef UTOPIA_TPL_ELASTICITY_GRADIENT_HPP
#define UTOPIA_TPL_ELASTICITY_GRADIENT_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

template <typename T>
UTOPIA_FUNCTION void elastic_material_gradient(const T mu,
                                               const T lmbda,
                                               const T *UTOPIA_RESTRICT f,
                                               const T *grad_test,
                                               const T dx,
                                               T *UTOPIA_RESTRICT lf,
                                               const int offset_i = 0)
{
        using namespace utopia::device;
        // Automatically generated
        T x0 = 22*lmbda*(11*lmbda*(-f[0]*f[3] + f[1]*f[2] + 1) + 6*mu);
T x1 = 2*f[0];
T x2 = 1.0/(pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + 1);
T x3 = 88*lmbda*mu;
T x4 = (1.0/132.0)/lmbda;
T x5 = grad_test[0]*x4;
T x6 = 2*f[1];
T x7 = grad_test[1]*x4;
T x8 = 2*f[2];
T x9 = 2*f[3];
lf[offset_i+0] += dx*(x5*(-f[3]*x0 + x3*(-x1*x2 + x1)) + x7*(f[2]*x0 + x3*(-x2*x6 + x6)));
lf[offset_i+1] += dx*(x5*(f[1]*x0 + x3*(-x2*x8 + x8)) + x7*(-f[0]*x0 + x3*(-x2*x9 + x9)));
}

#undef UTOPIA_RESTRICT

#endif  // UTOPIA_TPL_ELASTICITY_GRADIENT_HPP
