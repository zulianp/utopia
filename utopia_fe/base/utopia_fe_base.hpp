#ifndef UTOPIA_FE_BASE_HPP
#define UTOPIA_FE_BASE_HPP

#include "utopia_fe_config.hpp"
#include "utopia.hpp"

namespace utopia {
#ifdef WITH_TRILINOS_ALGEBRA
    using USparseMatrix = TpetraMatrix;
    using UVector  = TpetraVector;
#else
    using USparseMatrix = PetscMatrix;
    using UVector  = PetscVector;
#endif //WITH_TRILINOS_ALGEBRA

    using UIndexArray  = utopia::Traits<UVector>::IndexArray;
    using UScalarArray = utopia::Traits<UVector>::ScalarArray;
    using UIndexSet    = utopia::Traits<UVector>::IndexSet;

    using USerialMatrix = utopia::BlasMatrix<double>;
    using USerialVector = utopia::BlasVector<double>;
}

#endif //UTOPIA_FE_BASE_HPP
