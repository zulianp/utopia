#ifndef UTOPIA_INTREPID2_FORWARD_DECLARATIONS_HPP
#define UTOPIA_INTREPID2_FORWARD_DECLARATIONS_HPP

#include "utopia_kokkos_Field.hpp"

namespace utopia {

    namespace intrepid2 {

        // template <class Operator, typename Scalar>
        // class Assemble;

        template <typename Scalar>
        class FE;

        // template <typename Scalar>
        // class SubdomainValue;

        // template <typename FE>
        // class Field;

        // template <typename Scalar>
        // class FEAssembler;

    }  // namespace intrepid2

    template <typename T>
    class Traits<intrepid2::FE<T>> {
    public:
        using Scalar = T;
        using Field = utopia::kokkos::Field<intrepid2::FE<T>>;
    };
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FORWARD_DECLARATIONS_HPP