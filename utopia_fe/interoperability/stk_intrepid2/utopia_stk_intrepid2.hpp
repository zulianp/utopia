#ifndef UTOPIA_STK_INTREPID2_HPP
#define UTOPIA_STK_INTREPID2_HPP

#include "utopia_CreateFE.hpp"

#include "utopia_intrepid2_ForwardDeclarations.hpp"
#include "utopia_stk_ForwardDeclarations.hpp"

namespace utopia {

    template <typename Scalar>
    class CreateFE<utopia::stk::FunctionSpace, utopia::intrepid2::FE<Scalar>> {
    public:
        static void apply(const utopia::stk::FunctionSpace &space,
                          utopia::intrepid2::FE<Scalar> &fe,
                          const int degree = 0);
    };

}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_HPP