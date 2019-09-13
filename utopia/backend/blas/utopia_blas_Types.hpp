#ifndef UTOPIA_BLAS_TYPES_HPP
#define UTOPIA_BLAS_TYPES_HPP

#include "utopia_blas_Matrix.hpp"
#include "utopia_blas_Vector.hpp"
#include "utopia_blas_Traits.hpp"
#include "utopia_Tensor.hpp"

namespace utopia {
    using Matrixd = utopia::BlasDenseMatrix<double>;
    using Vectord = utopia::BlasVector<double>;
}

#endif //UTOPIA_BLAS_TYPES_HPP
