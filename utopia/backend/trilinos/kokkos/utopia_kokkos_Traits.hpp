#ifndef UTOPIA_KOKKOS_TRAITS_HPP
#define UTOPIA_KOKKOS_TRAITS_HPP

#include "utopia_kokkos_ForwardDeclarations.hpp"
#include "utopia_trilinos_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_BackendInfo.hpp"

namespace utopia {

    template<class Scalar_>
    class KokkosTraits {
    public:

        using Scalar   = Scalar_;
        using SizeType = utopia::TpetraSizeType;
        // using Matrix   = utopia::TpetraMatrix;
        // using SparseMatrix = utopia::TpetraMatrix;

        using KokkosView1 = Kokkos::View<Scalar *>;
        using Vector      = utopia::VectorView<KokkosView1>;

        //FIXME use Kokkos compatible wrapper
        using IndexSet    = utopia::TpetraIndexSet;
        using IndexArray  = utopia::TpetraIndexArray;
        using ScalarArray = utopia::TpetraScalarArray;

        enum {
            Backend = KOKKOS
        };

        static BackendInfo &backend_info()
        {
            static BackendInfo instance_("kokkos");
            return instance_;
        }
    };

    template<typename Scalar>
    using DefaultVectorView = utopia::VectorView<Kokkos::View<Scalar *>>;

    UTOPIA_MAKE_TRAITS_TPL_1(DefaultVectorView, KokkosTraits, 1);
}

#endif //UTOPIA_KOKKOS_TRAITS_HPP
