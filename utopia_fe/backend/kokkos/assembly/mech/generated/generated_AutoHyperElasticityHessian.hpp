#ifndef UTOPIA_TPL_ELASTICITY_HESSIAN_HPP
#define UTOPIA_TPL_ELASTICITY_HESSIAN_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

#ifndef AUTO_HYPER_ELASTICITY_DIM
#define AUTO_HYPER_ELASTICITY_DIM 2
#endif //AUTO_HYPER_ELASTICITY_DIM

template <typename T>
UTOPIA_FUNCTION void elastic_material_hessian(const T mu,
                                              const T lmbda,
                                              const T *UTOPIA_RESTRICT f,
                                              const T *grad_test,
                                              const T *grad_trial,
                                              const T dx,
                                              T *UTOPIA_RESTRICT bf,
                                              const int offset_ij = 0)
{
        using namespace utopia::device;
        // Automatically generated
        T x0 = 16*mu;
T x1 = f[0]*f[1];
T x2 = pow(f[0], 2);
T x3 = pow(f[1], 2);
T x4 = pow(f[2], 2);
T x5 = pow(f[3], 2);
T x6 = pow(x2 + x3 + x4 + x5 + 1, 2);
T x7 = 11*lmbda;
T x8 = x6*x7;
T x9 = f[2]*x8;
T x10 = -f[3]*x9 + x0*x1;
T x11 = -x3;
T x12 = -x4;
T x13 = -x5;
T x14 = x12 + x13 + x6 - 1;
T x15 = 8*mu;
T x16 = (1.0/6.0)/x6;
T x17 = grad_trial[0]*x16;
T x18 = -x2;
T x19 = grad_trial[1]*x16;
T x20 = f[2]*x0;
T x21 = f[1]*f[3];
T x22 = grad_test[0]*(f[0]*x20 - x21*x8);
T x23 = f[1]*f[2];
T x24 = 6*mu;
T x25 = f[0]*f[3];
T x26 = x7*(x23 - x25 + 1);
T x27 = x0*x23 + x6*(x23*x7 + x24 + x26);
T x28 = -grad_test[1]*(f[0]*x9 - x0*x21);
T x29 = x0*x25 + x6*(-x24 + x25*x7 - x26);
T x30 = -f[3]*x20 + x1*x8;
T x31 = x11 + x18 + x6 - 1;
bf[offset_ij+0] += dx*(x17*(grad_test[0]*(x15*(x11 + x14 + x2) + x5*x8) + grad_test[1]*x10) + x19*(grad_test[0]*x10 + grad_test[1]*(x15*(x14 + x18 + x3) + x4*x8)));
bf[offset_ij+1] += dx*(x17*(grad_test[1]*x27 + x22) + x19*(grad_test[0]*x29 + x28));
bf[offset_ij+2] += dx*(x17*(grad_test[1]*x29 + x22) + x19*(grad_test[0]*x27 + x28));
bf[offset_ij+3] += dx*(x17*(grad_test[0]*(x15*(x13 + x31 + x4) + x3*x8) - grad_test[1]*x30) + x19*(-grad_test[0]*x30 + grad_test[1]*(x15*(x12 + x31 + x5) + x2*x8)));
}

#undef UTOPIA_RESTRICT

#endif  // UTOPIA_TPL_ELASTICITY_HESSIAN_HPP
