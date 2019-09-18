#ifndef UTOPIA_TRILINOSTRAITS_HPP
#define UTOPIA_TRILINOSTRAITS_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_trilinos_ForwardDeclaration.hpp"
#include "utopia_trilinos_Base.hpp"
#include "utopia_BackendInfo.hpp"

#include <vector>

namespace utopia {

    class TpetraTraits {
    public:
        using Scalar   = utopia::TpetraScalar;
        using SizeType = utopia::TpetraSizeType;
        using Matrix   = utopia::TpetraMatrix;
        using SparseMatrix = utopia::TpetraMatrix;
        using Vector   = utopia::TpetraVector;

        //FIXME use Kokkos compatible wrapper
        using IndexSet    = utopia::TpetraIndexSet;
        using IndexArray  = utopia::TpetraIndexArray;
        using ScalarArray = utopia::TpetraScalarArray;

        enum {
            Backend = TRILINOS
        };

        static BackendInfo &backend_info()
        {
            static BackendInfo instance_("tpetra");
            return instance_;
        }
    };

    UTOPIA_MAKE_TRAITS_SPARSE(TpetraMatrix, TpetraTraits);
    UTOPIA_MAKE_TRAITS(TpetraVector, TpetraTraits);

    UTOPIA_MAKE_PARALLEL_TRAITS(TpetraVector);
    UTOPIA_MAKE_PARALLEL_TRAITS(TpetraMatrix);
}

#endif //UTOPIA_TRILINOSTRAITS_HPP
