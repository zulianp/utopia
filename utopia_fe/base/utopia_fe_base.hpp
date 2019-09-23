#ifndef UTOPIA_FE_BASE_HPP
#define UTOPIA_FE_BASE_HPP

#include "utopia_fe_config.hpp"
#include "utopia.hpp"

namespace utopia {
#ifdef WITH_TRILINOS_ALGEBRA
    using USparseMatrix = TSMatrixd;
    using UVector  = TVectord;
#else
    using USparseMatrix = DSMatrixd;
    using UVector  = DVectord;
#endif //WITH_TRILINOS_ALGEBRA

    using UIndexArray  = utopia::Traits<UVector>::IndexArray;
    using UScalarArray = utopia::Traits<UVector>::ScalarArray;
    using UIndexSet    = utopia::Traits<UVector>::IndexSet;

    using USerialMatrix = utopia::BlasDenseMatrix<double>;
    using USerialVector = utopia::BlasVector<double>;
}

#endif //UTOPIA_FE_BASE_HPP
