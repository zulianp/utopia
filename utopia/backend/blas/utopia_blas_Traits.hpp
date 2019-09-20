

#ifndef UTOPIA_UTOPIA_BLASTRAITS_HPP
#define UTOPIA_UTOPIA_BLASTRAITS_HPP

#include "utopia_Traits.hpp"
#include "utopia_blas_ForwardDeclarations.hpp"
#include <vector>

namespace utopia {

    template<typename T>
    class BLASTraits {
    public:
        using Scalar      = T;
        using Matrix      = utopia::BlasDenseMatrix<T>;
        using Vector      = utopia::BlasVector<T>;
        using SizeType    = std::size_t;
        using IndexSet    = utopia::BlasIndexSet;
        using ScalarArray = utopia::BlasArray<T>;
        using IndexArray  = utopia::BlasArray<int>;

        enum {
            Backend = BLAS
        };

        static BackendInfo &backend_info()
        {
            static BackendInfo instance_("blas");
            return instance_;
        }
    };

    UTOPIA_MAKE_TRAITS_TPL_1(BlasVector, BLASTraits, 1);
    UTOPIA_MAKE_TRAITS_DENSE_TPL_1(BlasDenseMatrix, BLASTraits, 2);

}

#endif //UTOPIA_UTOPIA_BLASTRAITS_HPP