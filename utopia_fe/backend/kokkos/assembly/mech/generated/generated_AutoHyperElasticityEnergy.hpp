#ifndef UTOPIA_TPL_ELASTICITY_ENERGY_HPP
#define UTOPIA_TPL_ELASTICITY_ENERGY_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

#ifndef AUTO_HYPER_ELASTICITY_DIM
#define AUTO_HYPER_ELASTICITY_DIM 2
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
        T x0 = pow(f[0], 2) + pow(f[1], 2) + pow(f[2], 2) + pow(f[3], 2);
e += (1.0/132.0)*dx*(88*lmbda*mu*(x0 - log(x0 + 1) - 2) + pow(11*lmbda*(-f[0]*f[3] + f[1]*f[2] + 1) + 6*mu, 2))/lmbda;
}

#undef UTOPIA_RESTRICT

#endif //UTOPIA_TPL_ELASTICITY_ENERGY_HPP
