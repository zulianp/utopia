#ifndef UTOPIA_INTREPID2_FORWARD_DECLARATIONS_HPP
#define UTOPIA_INTREPID2_FORWARD_DECLARATIONS_HPP

#include "utopia_kokkos_Field.hpp"

namespace utopia {

    namespace intrepid2 {
        template <typename Scalar>
        class FE;
    }  // namespace intrepid2

    template <typename T>
    class Traits<intrepid2::FE<T>> {
    public:
        using Scalar = T;
        using Field = utopia::kokkos::Field<intrepid2::FE<T>>;
    };
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FORWARD_DECLARATIONS_HPP