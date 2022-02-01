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
{{
        using namespace utopia::device;
        // Automatically generated
        {code}
}}

#undef UTOPIA_RESTRICT

#endif //UTOPIA_TPL_ELASTICITY_ENERGY_HPP
