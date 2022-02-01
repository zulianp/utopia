#ifndef UTOPIA_TPL_ELASTICITY_HPP
#define UTOPIA_TPL_ELASTICITY_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

template <typename T>
UTOPIA_FUNCTION void elastic_material(const T mu,
                                      const T lmbda,
                                      const T *UTOPIA_RESTRICT f,
                                      const T *grad_test,
                                      const T *grad_trial,
                                      const T dx,
                                      T &e,
                                      T *UTOPIA_RESTRICT lf,
                                      T *UTOPIA_RESTRICT bf,
                                      const int offset_i,
                                      const int offset_ij)
{
        using namespace utopia::device;
        // Automatically generated
        T x0 = f[0]*f[3];
T x1 = f[1]*f[2];
T x2 = x0 - x1;
T x3 = log(x2);
T x4 = lmbda*x3;
T x5 = mu - x4;
T x6 = lmbda + x5;
T x7 = f[3]*x6;
T x8 = f[2]*x7;
T x9 = pow(f[3], 2);
T x10 = lmbda*x9;
T x11 = pow(x2, 2);
T x12 = mu*x11;
T x13 = 1.0/x11;
T x14 = grad_trial[0]*x13;
T x15 = pow(f[2], 2);
T x16 = grad_trial[1]*x13;
T x17 = f[1]*grad_test[0];
T x18 = -x17*x7;
T x19 = x1*x6 + x2*x5;
T x20 = f[0]*x6;
T x21 = grad_test[1]*x20;
T x22 = -f[2]*x21;
T x23 = x0*x6 + x2*(-mu + x4);
T x24 = pow(f[1], 2);
T x25 = pow(f[0], 2);
T x26 = f[0]*mu;
T x27 = 1.0/x2;
T x28 = f[3]*mu;
T x29 = x27*x4;
T x30 = f[1]*mu;
T x31 = f[2]*mu;
bf[offset_ij+0] += dx*(x14*(grad_test[0]*(mu*x9 - x10*x3 + x10 + x12) - grad_test[1]*x8) + x16*(-grad_test[0]*x8 + grad_test[1]*(lmbda*x15 + mu*x15 + x12 - x15*x4)));
bf[offset_ij+1] += dx*(x14*(grad_test[1]*x19 + x18) + x16*(grad_test[0]*x23 + x22));
bf[offset_ij+2] += dx*(x14*(grad_test[1]*x23 + x18) + x16*(grad_test[0]*x19 + x22));
bf[offset_ij+3] += dx*(x14*(-f[1]*x21 + grad_test[0]*(lmbda*x24 + mu*x24 + x12 - x24*x4)) + x16*(grad_test[1]*(lmbda*x25 + mu*x25 + x12 - x25*x4) - x17*x20));
lf[offset_i+0] += dx*(grad_test[0]*(f[3]*x29 + x26 - x27*x28) + grad_test[1]*(-f[2]*x29 + x27*x31 + x30));
lf[offset_i+1] += dx*(grad_test[0]*(-f[1]*x29 + x27*x30 + x31) + grad_test[1]*(f[0]*x29 - x26*x27 + x28));
e += dx*((1.0/2.0)*lmbda*pow(x3, 2) - mu*x3 + (1.0/2.0)*mu*(x15 + x24 + x25 + x9 - 2));
}

#undef UTOPIA_RESTRICT

#endif //UTOPIA_TPL_ELASTICITY_HPP
