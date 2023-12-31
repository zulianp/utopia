

#ifndef UTOPIA_UTOPIA_BLASTRAITS_HPP
#define UTOPIA_UTOPIA_BLASTRAITS_HPP

#include "utopia_BackendInfo.hpp"
#include "utopia_Communicator.hpp"
#include "utopia_Device.hpp"
#include "utopia_Layout.hpp"
#include "utopia_Traits.hpp"
#include "utopia_blas_ForwardDeclarations.hpp"

#include <vector>

namespace utopia {

    template <typename T>
    class BLASTraits {
    public:
        using Scalar = T;
        using Matrix = utopia::BlasMatrix<T>;
        using Vector = utopia::BlasVector<T>;
        using SizeType = std::size_t;
        using LocalSizeType = std::size_t;
        using IndexSet = utopia::BlasIndexSet;
        using ScalarArray = utopia::BlasArray<T>;
        using IndexArray = utopia::BlasArray<int>;
        using Communicator = utopia::SelfCommunicator;
        using Device = utopia::Device<BLAS>;
        using Layout = utopia::Layout<SelfCommunicator, 1, SizeType>;
        using MatrixLayout = utopia::Layout<SelfCommunicator, 2, SizeType>;

        enum { Backend = BLAS };

        static BackendInfo &backend_info() {
            static BackendInfo instance_("blas");
            return instance_;
        }

        // tells petsc to device local size automatically
        static constexpr SizeType decide() { return 0; }

        // tells petsc to compute global size automatically
        static constexpr SizeType determine() { return 0; }
    };

    UTOPIA_MAKE_TRAITS_TPL_1(BlasVector, BLASTraits, 1);
    UTOPIA_MAKE_TRAITS_DENSE_TPL_1(BlasMatrix, BLASTraits, 2);

}  // namespace utopia

#endif  // UTOPIA_UTOPIA_BLASTRAITS_HPP