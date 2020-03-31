
#ifndef UTOPIA_TRILINOSTRAITS_HPP
#define UTOPIA_TRILINOSTRAITS_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_trilinos_ForwardDeclarations.hpp"
#include "utopia_trilinos_Base.hpp"
#include "utopia_BackendInfo.hpp"
#include "utopia_Device.hpp"
#include "utopia_Layout.hpp"

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
        using Layout        = utopia::Layout<TrilinosCommunicator, 1, LocalSizeType, SizeType>;
        using MatrixLayout  = utopia::Layout<TrilinosCommunicator, 2, LocalSizeType, SizeType>;

        enum {
            Backend = TRILINOS
        };

        static BackendInfo &backend_info()
        {
            static BackendInfo instance_("trilinos");
            return instance_;
        }

        //tells tpetra to device local size automatically
        static constexpr std::size_t decide()
        {
            return static_cast<std::size_t>(INVALID_INDEX);
        }

        //tells tpetra to compute global size automatically
        static constexpr Tpetra::global_size_t determine()
        {
            return Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()
        }
    };

    UTOPIA_MAKE_TRAITS_SPARSE(TpetraMatrix, TpetraTraits, 2);
    UTOPIA_MAKE_TRAITS(TpetraVector, TpetraTraits, 1);

    UTOPIA_MAKE_PARALLEL_TRAITS(TpetraVector);
    UTOPIA_MAKE_PARALLEL_TRAITS(TpetraMatrix);
}

#endif //UTOPIA_TRILINOSTRAITS_HPP
