#ifndef UTOPIA_KOKKOS_TRAITS_HPP
#define UTOPIA_KOKKOS_TRAITS_HPP

#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template<typename T, Size_t... Args>
    class Traits< ArrayView<T, Args...> > {
    public:
        using Scalar = T;
        using SizeType = Size_t;

        static const int Backend = SERIAL_HOMEMADE;
        static const int FILL_TYPE = FillType::DENSE;
    };

    template<class View>
    class Traits< VectorView<View> > : public Traits<View> {};

    template<class View>
    class Traits< MatrixView<View> > : public Traits<View> {};
}

#endif //UTOPIA_KOKKOS_TRAITS_HPP
