#ifndef UTOPIA_SADDLE_POINT_SOLVE_HPP
#define UTOPIA_SADDLE_POINT_SOLVE_HPP

namespace utopia {

    template <class FunctionSpace>
    class SaddlePointSolve : public Configurable {
    public:
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using FEFunctionInterface_t = utopia::FEFunctionInterface<FunctionSpace>;
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;

        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;

        using NewtonBase_t = utopia::NewtonBase<Matrix_t, Vector_t>;
        using Newton_t = utopia::Newton<Matrix_t, Vector_t>;

        class SubSystem {
        public:
            Matrix_t hessian;
            Vector_t solution, g;
        };

        void read(Input &in) override {
            auto space_from = std::make_shared<FunctionSpace>();
            auto space_to = std::make_shared<FunctionSpace>();

            solver_.resize(2);
            for (auto &s : solver_) {
                s = std::make_shared<OmniLinearSolver_t>();
            }

            in.get("picard_iterations", picard_iterations);
            in.get("first_order", first_order);
            in.get("dumping", dumping);
            in.get("write_matlab", write_matlab);
            in.get("augment_diagonal", augment_diagonal);

            in.require("from", [&](Input &node) {
                node.require("space", *space_from);
                auto fun_from = std::make_shared<FEModelFunction_t>(space_from);
                node.require("function", *fun_from);
                functions_.push_back(fun_from);

                node.get("solver", *solver_[0]);
            });

            in.require("to", [&](Input &node) {
                node.require("space", *space_to);
                auto fun_to = std::make_shared<FEModelFunction_t>(space_to);
                node.require("function", *fun_to);
                functions_.push_back(fun_to);

                node.get("solver", *solver_[1]);
            });

            in.get("transfer", transfer_);

            if (!transfer_.init(space_from, space_to)) {
                Utopia::Abort("FETransfer::init failed!");
            }
        }

        bool run() {
            std::size_t n = functions_.size();
            std::vector<SubSystem> systems(n);

            SubSystem condensed_system;

            bool hessians_constant = true;
            for (std::size_t i = 0; i < n; ++i) {
                auto &s = systems[i];
                auto &f = functions_[i];

                f->create_solution_vector(s.solution);
                f->hessian_and_gradient(s.solution, s.hessian, s.g);
                f->space()->apply_zero_constraints(s.g);

                rename("h_" + std::to_string(i), s.hessian);
                rename("g_" + std::to_string(i), s.g);

                if (i == 1) {
                    // Remove Dirichlet contributions from matrix
                    f->space()->apply_constraints(s.hessian, 0.);
                }

                if (write_matlab) {
                    write("H_" + std::to_string(i) + ".m", s.hessian);
                    write("G_" + std::to_string(i) + ".m", s.g);
                }

                if (!f->is_hessian_constant()) {
                    hessians_constant = false;
                }
            }

            rename("t", *transfer_.transfer_matrix());

            if (write_matlab) {
                write("T.m", *transfer_.transfer_matrix());
            }

            if (augment_diagonal && hessians_constant) {
                Matrix_t hessian_ptap;
                transfer_.apply(systems[1].hessian, hessian_ptap);

                //////////////////////////////////////////////////////////////////////
                condensed_system.hessian = systems[0].hessian + hessian_ptap;

                Vector_t gradient_p;
                transfer_.apply(systems[1].g, gradient_p);
                condensed_system.g = systems[0].g + gradient_p;
                functions_[0]->space()->apply_zero_constraints(condensed_system.g);
                //////////////////////////////////////////////////////////////////////

                Vector_t d = diag(hessian_ptap);
                systems[0].hessian += diag(d);
                functions_[0]->space()->comm().root_print("Diagonal 1 added to system 0!");
            }

            Vector_t g_coupled, diff_solution;
            Vector_t delta_from(layout(systems[0].g), 0.);
            Vector_t delta_to(layout(systems[1].g), 0.);
            Vector_t transferred(layout(systems[1].g), 0.);
            Vector_t lagrange_multiplier(systems[1].g);

            ////////////////////////////////////////////////////

            Scalar_t prev_norm = 1e6;

            IO<FunctionSpace> io(*functions_[0]->space());
            io.set_output_path("iterations.e");

            for (int k = 0; k < picard_iterations; ++k) {
                io.write(systems[0].solution, k, k);

                // Transfer response
                transfer_.apply_transpose(lagrange_multiplier, g_coupled);

                if (functions_[0]->is_hessian_constant()) {
                    systems[0].g.set(0.);
                    functions_[0]->gradient(systems[0].solution, systems[0].g);
                } else {
                    functions_[0]->hessian_and_gradient(systems[0].solution, systems[0].hessian, systems[0].g);
                }

                Scalar_t norm_lagrange_multiplier = norm2(g_coupled);

                // Add gradient contribution
                g_coupled += systems[0].g;

                // Negative gradient direction
                g_coupled = -g_coupled;

                if (g_coupled.has_nan_or_inf()) {
                    Utopia::Abort("g_coupled has NaN");
                }

                Scalar_t norm_g = norm2(g_coupled);

                if (norm_g > prev_norm) {
                    functions_[0]->space()->comm().root_print("Diverging!");
                    // break;
                }

                prev_norm = norm_g;

                functions_[0]->space()->comm().root_print("norm_g " + std::to_string(k) + ": " +
                                                          std::to_string(norm_g));

                functions_[0]->space()->comm().root_print("norm_Lagr " + std::to_string(k) + ": " +
                                                          std::to_string(norm_lagrange_multiplier));

                functions_[0]->space()->apply_zero_constraints(g_coupled);

                delta_from.set(0.);

                if (functions_[0]->is_hessian_constant() && k != 0) {
                    // Save cost of initialization of operator
                    solver_[0]->apply(g_coupled, delta_from);
                } else {
                    solver_[0]->solve(systems[0].hessian, g_coupled, delta_from);
                }

                if (delta_from.has_nan_or_inf()) {
                    Utopia::Abort("delta_from has NaN");
                }

                systems[0].solution += dumping * delta_from;

                //////////////////////////////////////////////////////////

                if (first_order) {
                    transfer_.apply(systems[0].solution, systems[1].solution);
                    // functions_[1]->space()->apply_constraints(systems[1].solution);

                    if (systems[1].solution.has_nan_or_inf()) {
                        Utopia::Abort("solution 1 has NaN");
                    }

                    ////////////////////////////////////

                    systems[1].g.set(0.);
                    functions_[1]->gradient(systems[1].solution, systems[1].g);
                } else {
                    transfer_.apply(systems[0].solution, diff_solution);

                    diff_solution -= systems[1].solution;

                    Scalar_t norm_diff = norm2(diff_solution);
                    functions_[0]->space()->comm().root_print("norm_diff: " + std::to_string(norm_diff));

                    transferred = systems[1].hessian * diff_solution;

                    functions_[1]->space()->apply_zero_constraints(transferred);

                    delta_to.set(0.);

                    if (functions_[1]->is_hessian_constant() && k != 0) {
                        solver_[1]->apply(transferred, delta_to);
                    } else {
                        solver_[1]->solve(systems[1].hessian, transferred, delta_to);
                    }

                    systems[1].solution += delta_to;

                    if (functions_[1]->is_hessian_constant()) {
                        functions_[1]->gradient(systems[1].solution, systems[1].g);
                    } else {
                        functions_[1]->hessian_and_gradient(systems[1].solution, systems[1].hessian, systems[1].g);
                    }
                }

                ////////////////////////////////////////////////////////

                // Clean boundary conditions
                // functions_[1]->space()->apply_zero_constraints(systems[1].g);

                // Copy g to lagrange multiplier
                lagrange_multiplier = systems[1].g;
            }

            ////////////////////////////////////////////////////////

            functions_[0]->report_solution(systems[0].solution);
            functions_[1]->report_solution(systems[1].solution);
            return true;
        }

        bool valid() const { return functions_.size() == 2; }

    private:
        std::vector<std::shared_ptr<FEFunctionInterface_t>> functions_;
        FETransfer<FunctionSpace> transfer_;
        std::vector<std::shared_ptr<OmniLinearSolver_t>> solver_;

        int picard_iterations{2};
        bool first_order{false};
        Scalar_t dumping{1};
        bool write_matlab{false};
        bool augment_diagonal{true};
    };
}  // namespace utopia

#endif  // UTOPIA_SADDLE_POINT_SOLVE_HPP
