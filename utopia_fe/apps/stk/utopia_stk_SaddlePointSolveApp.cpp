#include "utopia_Main.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_SaddlePointQPSolve.hpp"
#include "utopia_SaddlePointSolve.hpp"
#include "utopia_SemiGeometricMultigridNew.hpp"

#include "utopia_moonolith_stk_FETransfer.hpp"
#include "utopia_stk.hpp"

#include "utopia_stk_intrepid2.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

#include "utopia_moonolith_stk_Obstacle.hpp"

namespace utopia {
    template class NewmarkIntegrator<utopia::stk::FunctionSpace>;
    template class ImplicitEulerIntegrator<utopia::stk::FunctionSpace>;
}  // namespace utopia

void stk_saddle_point_solve(utopia::Input &in) {
    utopia::SaddlePointSolve<utopia::stk::FunctionSpace> app;
    app.read(in);

    if (app.valid()) {
        app.run();
    } else {
        utopia::err() << "stk_saddle_point_solve: invalid app setup\n";
    }
}

UTOPIA_REGISTER_APP(stk_saddle_point_solve);

void stk_saddle_point_qp_solve(utopia::Input &in) {
    utopia::SaddlePointQPSolve<utopia::stk::FunctionSpace> app;
    app.read(in);

    if (app.valid()) {
        app.run();
    } else {
        utopia::err() << "stk_saddle_point_qp_solve: invalid app setup\n";
    }
}

UTOPIA_REGISTER_APP(stk_saddle_point_qp_solve);
