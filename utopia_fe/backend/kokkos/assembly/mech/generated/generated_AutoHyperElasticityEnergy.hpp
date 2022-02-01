#ifndef UTOPIA_TPL_ELASTICITY_ENERGY_HPP
#define UTOPIA_TPL_ELASTICITY_ENERGY_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

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
        T x0 = log(f[0]*f[3] - f[1]*f[2]);
e += dx*((1.0/2.0)*lmbda*pow(x0, 2) - mu*x0 + (1.0/2.0)*mu*(pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2) - 2));
}

#undef UTOPIA_RESTRICT

#endif //UTOPIA_TPL_ELASTICITY_ENERGY_HPP
