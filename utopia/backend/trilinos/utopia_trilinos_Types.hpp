
#ifndef UTOPIA_TRILINOS_TYPES_HPP
#define UTOPIA_TRILINOS_TYPES_HPP

#include "utopia_Traits.hpp"
#include "utopia_trilinos_Traits.hpp"
#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_Vector.hpp"

namespace utopia {
    /*!
     * @brief Matrix representation of the trilinos backend.
     */

    using TSMatrixd = utopia::TpetraMatrix;
    using TVectord  = utopia::TpetraVector;
     
}

#endif //UTOPIA_UTOPIA_TRILINOS_TYPES_HPP

