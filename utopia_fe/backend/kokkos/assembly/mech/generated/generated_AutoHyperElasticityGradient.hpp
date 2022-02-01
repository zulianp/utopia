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
                                               const int offset_i)
{
        using namespace utopia::device;
        // Automatically generated
        T x0 = f[0]*mu;
T x1 = f[0]*f[3] - f[1]*f[2];
T x2 = 1.0/x1;
T x3 = f[3]*mu;
T x4 = lmbda*x2*log(x1);
T x5 = f[1]*mu;
T x6 = f[2]*mu;
lf[offset_i+0] += dx*(grad_test[0]*(f[3]*x4 + x0 - x2*x3) + grad_test[1]*(-f[2]*x4 + x2*x6 + x5));
lf[offset_i+1] += dx*(grad_test[0]*(-f[1]*x4 + x2*x5 + x6) + grad_test[1]*(f[0]*x4 - x0*x2 + x3));
}

#undef UTOPIA_RESTRICT

#endif  // UTOPIA_TPL_ELASTICITY_GRADIENT_HPP
