#ifndef UTOPIA_MORTAR_APP_HPP
#define UTOPIA_MORTAR_APP_HPP

#include "utopia_Main.hpp"

#include "utopia.hpp"

#include "utopia_SimpleNewton.hpp"

#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NCFunctionSpace.hpp"
#include "utopia_NCOmniAssembler.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_SemiGeometricMultigridNew.hpp"

#include "utopia_NLSolve.hpp"

namespace utopia {

    template <class FunctionSpace>
    class MortarApp {
    public:
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using NCFunctionSpace_t = utopia::NCFunctionSpace<FunctionSpace>;
        using FEModelFunction_t = utopia::FEModelFunction<NCFunctionSpace_t>;

        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;
        using Newton_t = utopia::Newton<Matrix_t, Vector_t>;

        void run(Input &in) {
            Communicator_t comm;
            NCFunctionSpace_t space(comm);
            in.get("space", space);

            auto problem = std::make_shared<FEModelFunction_t>(make_ref(space));
            in.get("problem", *problem);

            auto solver = std::make_shared<Newton_t>();
            solver->verbose(true);
            in.get("solver", *solver);

            NLSolve<NCFunctionSpace_t> nlsolve;
            nlsolve.init(problem);
            nlsolve.set_solver(solver);
            nlsolve.solve();
        }
    };
}  // namespace utopia

#endif  // UTOPIA_MORTAR_APP_HPP
