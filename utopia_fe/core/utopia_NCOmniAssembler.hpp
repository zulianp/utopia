#ifndef UTOPIA_NC_OMNI_ASSEMBLER_HPP
#define UTOPIA_NC_OMNI_ASSEMBLER_HPP

#include "utopia_FECoreForwardDeclarations.hpp"
#include "utopia_NCFunctionSpace.hpp"

namespace utopia {

    template <class FunctionSpace>
    class OmniAssembler<NCFunctionSpace<FunctionSpace>> final : public OmniAssembler<FunctionSpace> {
    public:
        using Super = utopia::OmniAssembler<FunctionSpace>;
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using Environment = utopia::Environment<NCFunctionSpace<FunctionSpace>>;

        ~OmniAssembler() {}
        OmniAssembler(const std::shared_ptr<NCFunctionSpace<FunctionSpace>> &space)
            : Super(space->unconstrained_space()) {}

        void set_environment(const std::shared_ptr<Environment> &) {
            utopia::err()
                << "OmniAssembler<NCFunctionSpace<FunctionSpace>>[Warning] set_environment not implemented!\n";
        }
    };

}  // namespace utopia

#endif  // UTOPIA_NC_OMNI_ASSEMBLER_HPP
