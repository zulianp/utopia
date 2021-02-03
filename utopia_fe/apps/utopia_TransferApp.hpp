#ifndef UTOPIA_TRANSFER_APP
#define UTOPIA_TRANSFER_APP

#include "utopia_FEApp.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_libmesh.hpp"

#include <memory>
#include <string>

namespace utopia {
    class LocalAssembler;
    class Local2Global;

    class TransferApp final : public FEApp {
    public:
        class InputSpace;

        ~TransferApp();
        TransferApp();

        void run(Input &in) override;

        static std::string command() { return "-transfer"; }

    private:
#ifdef UTOPIA_WITH_TINY_EXPR
        std::shared_ptr<SymbolicFunction> fun;
#else
        std::shared_ptr<ConstantCoefficient<double, 0> > fun;
#endif  // UTOPIA_WITH_TINY_EXPR

        bool fun_is_constant;
        bool write_operators_to_disk;
    };
}  // namespace utopia

#endif  // UTOPIA_TRANSFER_APP
