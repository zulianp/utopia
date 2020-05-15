#ifndef UTOPIA_BLAS_TYPES_HPP
#define UTOPIA_BLAS_TYPES_HPP

#include "utopia_Tensor.hpp"
#include "utopia_blas_Matrix.hpp"
#include "utopia_blas_Traits.hpp"
#include "utopia_blas_Vector.hpp"

namespace utopia {
    using BlasMatrixd = utopia::BlasMatrix<double>;
    using BlasVectord = utopia::BlasVector<double>;
}  // namespace utopia

#endif  // UTOPIA_BLAS_TYPES_HPP
