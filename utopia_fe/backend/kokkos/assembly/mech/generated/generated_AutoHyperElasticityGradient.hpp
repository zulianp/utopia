#ifndef UTOPIA_TPL_ELASTICITY_GRADIENT_HPP
#define UTOPIA_TPL_ELASTICITY_GRADIENT_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

#ifndef AUTO_HYPER_ELASTICITY_DIM
#define AUTO_HYPER_ELASTICITY_DIM 3
#endif //AUTO_HYPER_ELASTICITY_DIM

template <typename T>
UTOPIA_FUNCTION void elastic_material_gradient(const T mu,
                                               const T lmbda,
                                               const T *UTOPIA_RESTRICT f,
                                               const T *grad_test,
                                               const T dx,
                                               T *UTOPIA_RESTRICT lf)
{
        using namespace utopia::device;
        // Automatically generated
        T x0 = f[4]*f[8];
T x1 = f[5]*f[7];
T x2 = x0 - x1;
T x3 = f[5]*f[6];
T x4 = f[3]*f[7];
T x5 = f[3]*f[8];
T x6 = f[4]*f[6];
T x7 = f[0]*x0 - f[0]*x1 + f[1]*x3 - f[1]*x5 + f[2]*x4 - f[2]*x6;
T x8 = mu*x7;
T x9 = lmbda*log(x7);
T x10 = -x3 + x5;
T x11 = x4 - x6;
T x12 = dx/x7;
T x13 = f[1]*f[8] - f[2]*f[7];
T x14 = f[0]*f[8] - f[2]*f[6];
T x15 = f[0]*f[7] - f[1]*f[6];
T x16 = f[1]*f[5] - f[2]*f[4];
T x17 = f[0]*f[5] - f[2]*f[3];
T x18 = f[0]*f[4] - f[1]*f[3];
lf[0] += x12*(grad_test[0]*(f[0]*x8 - mu*x2 + x2*x9) + grad_test[1]*(f[1]*x8 + mu*x10 - x10*x9) + grad_test[2]*(f[2]*x8 - mu*x11 + x11*x9));
lf[1] += x12*(grad_test[0]*(f[3]*x8 + mu*x13 - x13*x9) + grad_test[1]*(f[4]*x8 - mu*x14 + x14*x9) + grad_test[2]*(f[5]*x8 + mu*x15 - x15*x9));
lf[2] += x12*(grad_test[0]*(f[6]*x8 - mu*x16 + x16*x9) + grad_test[1]*(f[7]*x8 + mu*x17 - x17*x9) + grad_test[2]*(f[8]*x8 - mu*x18 + x18*x9));
}

#undef UTOPIA_RESTRICT

#endif  // UTOPIA_TPL_ELASTICITY_GRADIENT_HPP
