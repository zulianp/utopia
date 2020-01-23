#ifndef UTOPIA_TRILINOSTRAITS_HPP
#define UTOPIA_TRILINOSTRAITS_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_trilinos_ForwardDeclarations.hpp"
#include "utopia_trilinos_Base.hpp"
#include "utopia_BackendInfo.hpp"
#include "utopia_Device.hpp"

#include <vector>

namespace utopia {

    class TpetraTraits {
    public:
        using Scalar        = utopia::TpetraScalar;
        using SizeType      = utopia::TpetraSizeType;
        using LocalSizeType = utopia::TpetraLocalSizeType;
        using Matrix        = utopia::TpetraMatrix;
        using SparseMatrix  = utopia::TpetraMatrix;
        using Vector        = utopia::TpetraVector;

        using IndexSet      = utopia::TpetraIndexSet;
        using IndexArray    = utopia::TpetraIndexArray;
        using ScalarArray   = utopia::TpetraScalarArray;

        using Communicator  = utopia::TrilinosCommunicator;
        using Node          = utopia::DefaultKokkosNode;
        using Device        = utopia::Device<TRILINOS>;

        enum {
            Backend = TRILINOS
        };

        static BackendInfo &backend_info()
        {
            static BackendInfo instance_("trilinos");
            return instance_;
        }
    };

    UTOPIA_MAKE_TRAITS_SPARSE(TpetraMatrix, TpetraTraits, 2);
    UTOPIA_MAKE_TRAITS(TpetraVector, TpetraTraits, 1);

    UTOPIA_MAKE_PARALLEL_TRAITS(TpetraVector);
    UTOPIA_MAKE_PARALLEL_TRAITS(TpetraMatrix);
}

#endif //UTOPIA_TRILINOSTRAITS_HPP
