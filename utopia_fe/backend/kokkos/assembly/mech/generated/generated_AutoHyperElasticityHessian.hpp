#ifndef UTOPIA_TPL_ELASTICITY_HESSIAN_HPP
#define UTOPIA_TPL_ELASTICITY_HESSIAN_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

#ifndef AUTO_HYPER_ELASTICITY_DIM
#define AUTO_HYPER_ELASTICITY_DIM 3
#endif //AUTO_HYPER_ELASTICITY_DIM

template <typename T>
UTOPIA_FUNCTION void elastic_material_hessian(const T mu,
                                              const T lmbda,
                                              const T *UTOPIA_RESTRICT f,
                                              const T *grad_test,
                                              const T *grad_trial,
                                              const T dx,
                                              T *UTOPIA_RESTRICT bf)
{
        using namespace utopia::device;
        // Automatically generated
        T x0 = f[4]*f[8];
T x1 = f[5]*f[6];
T x2 = f[3]*f[7];
T x3 = f[5]*f[7];
T x4 = f[3]*f[8];
T x5 = f[4]*f[6];
T x6 = f[0]*x0 - f[0]*x3 + f[1]*x1 - f[1]*x4 + f[2]*x2 - f[2]*x5;
T x7 = log(x6);
T x8 = lmbda*x7;
T x9 = mu - x8;
T x10 = lmbda + x9;
T x11 = x0 - x3;
T x12 = x10*x11;
T x13 = x2 - x5;
T x14 = grad_test[2]*x13;
T x15 = -x1 + x4;
T x16 = x12*x15;
T x17 = pow(x11, 2);
T x18 = lmbda*x17;
T x19 = pow(x6, 2);
T x20 = mu*x19;
T x21 = x10*x15;
T x22 = pow(x13, 2);
T x23 = pow(x15, 2);
T x24 = dx/x19;
T x25 = f[1]*f[8] - f[2]*f[7];
T x26 = x10*x25;
T x27 = grad_test[0]*x11;
T x28 = -x26*x27;
T x29 = x6*x9;
T x30 = f[8]*x29 + x15*x26;
T x31 = f[7]*x29 + x13*x26;
T x32 = f[0]*f[8] - f[2]*f[6];
T x33 = x10*x32;
T x34 = grad_test[1]*x15;
T x35 = -x33*x34;
T x36 = f[6]*x29 + x13*x33;
T x37 = x6*(-mu + x8);
T x38 = f[8]*x37 + x11*x33;
T x39 = f[0]*f[7] - f[1]*f[6];
T x40 = x10*x39;
T x41 = -x14*x40;
T x42 = f[6]*x37 + x15*x40;
T x43 = f[7]*x37 + x11*x40;
T x44 = f[1]*f[5] - f[2]*f[4];
T x45 = x10*x44;
T x46 = x27*x45;
T x47 = f[4]*x29 + x13*x45;
T x48 = f[5]*x29 + x15*x45;
T x49 = f[0]*f[5] - f[2]*f[3];
T x50 = x10*x49;
T x51 = x34*x50;
T x52 = f[3]*x29 + x13*x50;
T x53 = f[5]*x37 + x11*x50;
T x54 = f[0]*f[4] - f[1]*f[3];
T x55 = x10*x54;
T x56 = x14*x55;
T x57 = f[4]*x37 + x11*x55;
T x58 = f[3]*x37 + x15*x55;
T x59 = grad_test[2]*x39;
T x60 = x26*x32;
T x61 = pow(x25, 2);
T x62 = grad_test[0]*x26;
T x63 = grad_test[1]*x33;
T x64 = pow(x39, 2);
T x65 = pow(x32, 2);
T x66 = -x44*x62;
T x67 = f[2]*x29 + x33*x44;
T x68 = f[1]*x29 + x40*x44;
T x69 = -x49*x63;
T x70 = f[0]*x29 + x40*x49;
T x71 = f[2]*x37 + x26*x49;
T x72 = grad_test[2]*x54;
T x73 = -x40*x72;
T x74 = f[0]*x37 + x33*x54;
T x75 = f[1]*x37 + x26*x54;
T x76 = x45*x49;
T x77 = pow(x44, 2);
T x78 = pow(x54, 2);
T x79 = pow(x49, 2);
bf[0] += x24*(grad_trial[0]*(grad_test[0]*(mu*x17 - x18*x7 + x18 + x20) - grad_test[1]*x16 + x12*x14) + grad_trial[1]*(-grad_test[0]*x16 + grad_test[1]*(lmbda*x23 + mu*x23 + x20 - x23*x8) - x14*x21) + grad_trial[2]*(grad_test[0]*x12*x13 - grad_test[1]*x13*x21 + grad_test[2]*(lmbda*x22 + mu*x22 + x20 - x22*x8)));
bf[1] += x24*(grad_trial[0]*(grad_test[1]*x30 - grad_test[2]*x31 + x28) + grad_trial[1]*(grad_test[0]*x38 + grad_test[2]*x36 + x35) + grad_trial[2]*(-grad_test[0]*x43 + grad_test[1]*x42 + x41));
bf[2] += x24*(grad_trial[0]*(-grad_test[1]*x48 + grad_test[2]*x47 + x46) + grad_trial[1]*(-grad_test[0]*x53 - grad_test[2]*x52 + x51) + grad_trial[2]*(grad_test[0]*x57 - grad_test[1]*x58 + x56));
bf[3] += x24*(grad_trial[0]*(grad_test[1]*x38 - grad_test[2]*x43 + x28) + grad_trial[1]*(grad_test[0]*x30 + grad_test[2]*x42 + x35) + grad_trial[2]*(-grad_test[0]*x31 + grad_test[1]*x36 + x41));
bf[4] += x24*(grad_trial[0]*(grad_test[0]*(lmbda*x61 + mu*x61 + x20 - x61*x8) - grad_test[1]*x60 + x26*x59) + grad_trial[1]*(-grad_test[0]*x60 + grad_test[1]*(lmbda*x65 + mu*x65 + x20 - x65*x8) - x33*x59) + grad_trial[2]*(grad_test[2]*(lmbda*x64 + mu*x64 + x20 - x64*x8) + x39*x62 - x39*x63));
bf[5] += x24*(grad_trial[0]*(grad_test[1]*x67 - grad_test[2]*x68 + x66) + grad_trial[1]*(grad_test[0]*x71 + grad_test[2]*x70 + x69) + grad_trial[2]*(-grad_test[0]*x75 + grad_test[1]*x74 + x73));
bf[6] += x24*(grad_trial[0]*(-grad_test[1]*x53 + grad_test[2]*x57 + x46) + grad_trial[1]*(-grad_test[0]*x48 - grad_test[2]*x58 + x51) + grad_trial[2]*(grad_test[0]*x47 - grad_test[1]*x52 + x56));
bf[7] += x24*(grad_trial[0]*(grad_test[1]*x71 - grad_test[2]*x75 + x66) + grad_trial[1]*(grad_test[0]*x67 + grad_test[2]*x74 + x69) + grad_trial[2]*(-grad_test[0]*x68 + grad_test[1]*x70 + x73));
bf[8] += x24*(grad_trial[0]*(grad_test[0]*(lmbda*x77 + mu*x77 + x20 - x77*x8) - grad_test[1]*x76 + x45*x72) + grad_trial[1]*(-grad_test[0]*x76 + grad_test[1]*(lmbda*x79 + mu*x79 + x20 - x79*x8) - x50*x72) + grad_trial[2]*(grad_test[0]*x45*x54 - grad_test[1]*x50*x54 + grad_test[2]*(lmbda*x78 + mu*x78 + x20 - x78*x8)));
}

#undef UTOPIA_RESTRICT

#endif  // UTOPIA_TPL_ELASTICITY_HESSIAN_HPP
