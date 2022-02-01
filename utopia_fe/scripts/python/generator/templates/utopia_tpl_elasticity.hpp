#ifndef UTOPIA_TPL_ELASTICITY_HPP
#define UTOPIA_TPL_ELASTICITY_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

template <typename T>
UTOPIA_FUNCTION void elastic_material(const T mu,
                                      const T lmbda,
                                      const T *UTOPIA_RESTRICT f,
                                      const T *grad_i,
                                      const T *grad_j,
                                      const T dx,
                                      T &e,
                                      T *UTOPIA_RESTRICT lf,
                                      T *UTOPIA_RESTRICT bf,
                                      const int offset_i,
                                      const int offset_ij)
{{
        using namespace utopia::device;
        // Automatically generated
        {code}
}}

#undef UTOPIA_RESTRICT

#endif //UTOPIA_TPL_ELASTICITY_HPP
