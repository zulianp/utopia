#include "utopia_Main.hpp"
#include "utopia_Multiphysics.hpp"

#include "utopia_NewmarkIntegrator.hpp"

#include "utopia_stk.hpp"

#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_moonolith_stk_Obstacle.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

namespace utopia {

    template class NewmarkIntegrator<utopia::stk::FunctionSpace>;

    template <class FunctionSpace>
    class NLSolveApp : public Configurable {
    public:
        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;

        void read(Input &in) override {}
        void run() {}

        // NewmarkIntegrator<FunctionSpace> integrator_;
    };

}  // namespace utopia

void stk_nlsolve(utopia::Input &in) {
    utopia::NLSolveApp<utopia::stk::FunctionSpace> app;
    app.read(in);
    app.run();
}

UTOPIA_REGISTER_APP(stk_nlsolve);
