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
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
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
            int max_it = 10;
            Scalar_t stol = 1e-6;
            bool verbose = true;
            bool debug_nonlinear_material = false;
            OmniLinearSolver_t solver;
            in.get("solver", [&](Input &in) {
                in.get("max_it", max_it);
                in.get("stol", stol);
                in.get("verbose", verbose);
                in.get("debug_nonlinear_material", debug_nonlinear_material);
                in.get("linear_solver", solver);
            });

            Vector_t x, g, c;
            Matrix_t H;

            space.create_matrix(H);
            space.create_vector(x);
            space.create_vector(g);

            c.zeros(layout(x));

            // Poor's man Newton solver
            for (int it = 0; it < max_it; ++it) {
                if (!assembler->assemble(x, H, g)) {
                    utopia::err() << "Failed to assemble!\n";
                    return;
                }

                g *= -1.0;

                if (it == 0) {
                    space.apply_constraints(H, g);
                } else {
                    space.apply_constraints(H);
                    space.apply_zero_constraints(g);
                }

                c.set(0.0);
                solver.solve(H, g, c);

                x += c;

                if (assembler->is_linear()) {
                    break;
                }

                if (debug_nonlinear_material) {
                    space.write("it_" + std::to_string(it) + ".e", x);
                }

                Scalar_t norm_c = norm2(c);

                if (verbose) {
                    utopia::out() << "norm_c: " << norm_c << '\n';
                }

                if (norm_c < stol) {
                    break;
                }
            }

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
