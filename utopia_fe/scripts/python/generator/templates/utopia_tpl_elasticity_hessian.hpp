#ifndef UTOPIA_TPL_ELASTICITY_HESSIAN_HPP
#define UTOPIA_TPL_ELASTICITY_HESSIAN_HPP

#define UTOPIA_RESTRICT __restrict__

#include "utopia_Algorithms.hpp"

template <typename T>
UTOPIA_FUNCTION void elastic_material_hessian(const T mu,
                                              const T lmbda,
                                              const T *UTOPIA_RESTRICT f,
                                              const T *grad_test,
                                              const T *grad_trial,
                                              const T dx,
                                              T *UTOPIA_RESTRICT bf,
                                              const int offset_ij = 0)
{{
        using namespace utopia::device;
        // Automatically generated
        {code}
}}

#undef UTOPIA_RESTRICT

#endif  // UTOPIA_TPL_ELASTICITY_HESSIAN_HPP
