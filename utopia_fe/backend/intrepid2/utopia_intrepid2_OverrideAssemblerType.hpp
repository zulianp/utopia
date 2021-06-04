#ifndef UTOPIA_INTREPID2_OVERRIDE_ASSEMBLER_TYPE_HPP
#define UTOPIA_INTREPID2_OVERRIDE_ASSEMBLER_TYPE_HPP

#include "utopia_Traits.hpp"
#include "utopia_intrepid2_ForwardDeclarations.hpp"

namespace utopia {
    namespace intrepid2 {

        template <class FunctionSpace_, class Material_>
        class OverrideAssemblerType {
        public:
            using Material = Material_;
            using FunctionSpace = FunctionSpace_;

            using Scalar = typename Traits<FunctionSpace>::Scalar;
            using Type = utopia::intrepid2::Assemble<Material_, Scalar>;
        };

    }  // namespace intrepid2

}  // namespace utopia

#endif  // UTOPIA_INTREPID2_OVERRIDE_ASSEMBLER_TYPE_HPP
