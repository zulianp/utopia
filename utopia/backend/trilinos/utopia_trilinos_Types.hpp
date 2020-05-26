
#ifndef UTOPIA_TRILINOS_TYPES_HPP
#define UTOPIA_TRILINOS_TYPES_HPP

#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_Vector.hpp"
#include "utopia_Traits.hpp"
#include "utopia_trilinos_Traits.hpp"

namespace utopia {
    /*!
     * @brief Matrix representation of the trilinos backend.
     */

    using TpetraMatrixd = utopia::TpetraMatrix;
    using TpetraVectord = utopia::TpetraVector;
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_TRILINOS_TYPES_HPP
