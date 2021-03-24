#ifndef UTOPIA_STK_INTREPID2_OMNI_ASSEMBLER_HPP
#define UTOPIA_STK_INTREPID2_OMNI_ASSEMBLER_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_stk_ForwardDeclarations.hpp"
#include "utopia_stk_FunctionSpace.hpp"

namespace utopia {
    namespace stk {

        class OmniAssembler : public Configurable {
        public:
            using Matrix = Traits<stk::FunctionSpace>::Matrix;
            using Vector = Traits<stk::FunctionSpace>::Vector;
            using Scalar = Traits<stk::FunctionSpace>::Scalar;

            OmniAssembler(const std::shared_ptr<stk::FunctionSpace> &space);
            virtual ~OmniAssembler();

            bool assemble(const Vector &x, Matrix &jacobian, Vector &fun);
            void read(Input &in) override;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
        };
    }  // namespace stk

    template <>
    class OmniAssembler<utopia::stk::FunctionSpace> final : public utopia::stk::OmniAssembler {
    public:
        using Super = utopia::stk::OmniAssembler;
        using Super::Super;
    };

}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_OMNI_ASSEMBLER_HPP
