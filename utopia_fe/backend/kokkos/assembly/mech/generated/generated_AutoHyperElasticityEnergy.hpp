#ifndef UTOPIA_TPL_ELASTICITY_ENERGY_HPP
#define UTOPIA_TPL_ELASTICITY_ENERGY_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

#ifndef AUTO_HYPER_ELASTICITY_DIM
#define AUTO_HYPER_ELASTICITY_DIM 3
#endif //AUTO_HYPER_ELASTICITY_DIM

template <typename T>
UTOPIA_FUNCTION void elastic_material_energy(const T mu,
                                      const T lmbda,
                                      const T *UTOPIA_RESTRICT f,
                                      const T dx,
                                      T &e
                                    )
{
        using namespace utopia::device;
        // Automatically generated
        T x0 = log(f[0]*f[4]*f[8] - f[0]*f[5]*f[7] - f[1]*f[3]*f[8] + f[1]*f[5]*f[6] + f[2]*f[3]*f[7] - f[2]*f[4]*f[6]);
e += dx*((1.0/2.0)*lmbda*pow(x0, 2) - mu*x0 + (1.0/2.0)*mu*(pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) + pow(f[4], 2) + pow(f[5], 2) + pow(f[6], 2) + pow(f[7], 2) + pow(f[8], 2) - 3));
}

#undef UTOPIA_RESTRICT

#endif //UTOPIA_TPL_ELASTICITY_ENERGY_HPP
