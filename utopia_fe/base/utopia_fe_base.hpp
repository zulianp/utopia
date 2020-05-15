#ifndef UTOPIA_FE_BASE_HPP
#define UTOPIA_FE_BASE_HPP

#include "utopia.hpp"
#include "utopia_CoreDecprecatedHeaders.hpp"
#include "utopia_fe_config.hpp"

namespace utopia {
#ifdef WITH_TRILINOS_ALGEBRA
    using USparseMatrix = TpetraMatrixd;
    using UVector = TpetraVectord;
#else
    using USparseMatrix = PetscMatrix;
    using UVector = PetscVector;
#endif  // WITH_TRILINOS_ALGEBRA

    using UIndexArray = utopia::Traits<UVector>::IndexArray;
    using UScalarArray = utopia::Traits<UVector>::ScalarArray;
    using UIndexSet = utopia::Traits<UVector>::IndexSet;
    using UScalar = utopia::Traits<UVector>::Scalar;
    using USizeType = utopia::Traits<UVector>::SizeType;

    using USerialMatrix = utopia::BlasMatrix<UScalar>;
    using USerialVector = utopia::BlasVector<UScalar>;
}  // namespace utopia

#endif  // UTOPIA_FE_BASE_HPP
