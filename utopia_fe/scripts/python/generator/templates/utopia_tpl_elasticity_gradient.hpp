#ifndef UTOPIA_TPL_ELASTICITY_GRADIENT_HPP
#define UTOPIA_TPL_ELASTICITY_GRADIENT_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

template <typename T>
UTOPIA_FUNCTION void elastic_material_gradient(const T mu,
                                      const T lmbda,
                                      const T *UTOPIA_RESTRICT f,
                                      const T *grad_j,
                                      const T dx,
                                      T *UTOPIA_RESTRICT lf,
                                      const int offset_i
                                      )
{{
        using namespace utopia::device;
        // Automatically generated
        {code}
}}

#undef UTOPIA_RESTRICT

#endif //UTOPIA_TPL_ELASTICITY_GRADIENT_HPP
