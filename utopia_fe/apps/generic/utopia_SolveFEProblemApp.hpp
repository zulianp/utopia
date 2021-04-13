#ifndef UTOPIA_SOLVE_FE_PROBLEM_APP_HPP
#define UTOPIA_SOLVE_FE_PROBLEM_APP_HPP

#include "utopia_fe_base.hpp"

#include "utopia_Input.hpp"
#include "utopia_fe_Core.hpp"

namespace utopia {

    template <class FunctionSpace>
    class SolveFEProblemApp {
    public:
        // Extract front-end associated objects
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;

        // Use specialized compoenents for function space
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;

        void run(Input &in) {
            in.get("space", space);

            if (space.empty()) {
                utopia::err() << "[Error] space is undefined\n";
                return;
            }

            assembler = std::make_shared<OmniAssembler_t>(make_ref(space));
            in.get("assembly", *assembler);

            // nonlinear_solver = std::make_shared<Newton<Matrix_t, Vector_t>>(std::make_shared<OmniLinearSolver_t>());
            // in.get("solver", *nonlinear_solver);

            // REMOVE ME Once we have a nonlinear solver
            OmniLinearSolver_t solver;
            in.get("solver", [&](Input &in) { in.get("linear_solver", solver); });

            Vector_t x, g;
            Matrix_t H;

            space.create_matrix(H);
            space.create_vector(x);
            space.create_vector(g);
            if (!assembler->assemble(x, H, g)) {
                utopia::err() << "Failed to assemble!\n";
                return;
            }

            g *= -1.0;
            space.apply_constraints(H, g);

            solver.solve(H, g, x);

            std::string output = "output.e";
            in.get("output", output);
            space.write(output, x);
        }

        FunctionSpace space;
        std::shared_ptr<OmniAssembler_t> assembler;
        // std::shared_ptr<NonLinearSolver<Matrix_t, Vector_t>> nonlinear_solver;
    };
}  // namespace utopia

#endif  // UTOPIA_SOLVE_FE_PROBLEM_APP_HPP
