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
                                      const int offset_i = 0,
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
T x6 = x2 + x3 + x4 + x5;
T x7 = x6 + 1;
T x8 = pow(x7, 2);
T x9 = 11*lmbda;
T x10 = x8*x9;
T x11 = f[2]*x10;
T x12 = -f[3]*x11 + x0*x1;
T x13 = -x3;
T x14 = -x4;
T x15 = -x5;
T x16 = x14 + x15 + x8 - 1;
T x17 = 8*mu;
T x18 = (1.0/6.0)/x8;
T x19 = grad_trial[0]*x18;
T x20 = -x2;
T x21 = grad_trial[1]*x18;
T x22 = f[2]*x0;
T x23 = f[1]*f[3];
T x24 = grad_test[0]*(f[0]*x22 - x10*x23);
T x25 = f[1]*f[2];
T x26 = 6*mu;
T x27 = f[0]*f[3];
T x28 = x9*(x25 - x27 + 1);
T x29 = x26 + x28;
T x30 = x0*x25 + x8*(x25*x9 + x29);
T x31 = -grad_test[1]*(f[0]*x11 - x0*x23);
T x32 = x0*x27 + x8*(-x26 + x27*x9 - x28);
T x33 = -f[3]*x22 + x1*x10;
T x34 = x13 + x20 + x8 - 1;
T x35 = 22*lmbda*x29;
T x36 = 2*f[0];
T x37 = 1.0/x7;
T x38 = 88*lmbda*mu;
T x39 = (1.0/132.0)/lmbda;
T x40 = grad_test[0]*x39;
T x41 = 2*f[1];
T x42 = grad_test[1]*x39;
T x43 = 2*f[2];
T x44 = 2*f[3];
bf[offset_ij+0] += dx*(x19*(grad_test[0]*(x10*x5 + x17*(x13 + x16 + x2)) + grad_test[1]*x12) + x21*(grad_test[0]*x12 + grad_test[1]*(x10*x4 + x17*(x16 + x20 + x3))));
bf[offset_ij+1] += dx*(x19*(grad_test[1]*x30 + x24) + x21*(grad_test[0]*x32 + x31));
bf[offset_ij+2] += dx*(x19*(grad_test[1]*x32 + x24) + x21*(grad_test[0]*x30 + x31));
bf[offset_ij+3] += dx*(x19*(grad_test[0]*(x10*x3 + x17*(x15 + x34 + x4)) - grad_test[1]*x33) + x21*(-grad_test[0]*x33 + grad_test[1]*(x10*x2 + x17*(x14 + x34 + x5))));
lf[offset_i+0] += dx*(x40*(-f[3]*x35 + x38*(-x36*x37 + x36)) + x42*(f[2]*x35 + x38*(-x37*x41 + x41)));
lf[offset_i+1] += dx*(x40*(f[1]*x35 + x38*(-x37*x43 + x43)) + x42*(-f[0]*x35 + x38*(-x37*x44 + x44)));
e += dx*x39*(pow(x29, 2) + x38*(x6 - log(x7) - 2));
}

#undef UTOPIA_RESTRICT

#endif //UTOPIA_TPL_ELASTICITY_HPP
