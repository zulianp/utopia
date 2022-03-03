#ifndef UTOPIA_SADDLE_POINT_QP_SOLVE_HPP
#define UTOPIA_SADDLE_POINT_QP_SOLVE_HPP

#include "utopia_ActiveSet.hpp"
#include "utopia_polymorphic_QPSolver.hpp"

namespace utopia {

    template <class FunctionSpace>
    class SaddlePointQPSolve : public Configurable {
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
        using OmniQPSolver_t = utopia::OmniQPSolver<Matrix_t, Vector_t>;

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

            linear_solver_ = std::make_shared<OmniLinearSolver_t>();
            qp_solver_ = std::make_shared<OmniQPSolver_t>();

            in.get("picard_iterations", picard_iterations);
            in.get("first_order", first_order);
            in.get("dumping", dumping);
            in.get("write_matlab", write_matlab);
            in.get("augment_diagonal", augment_diagonal);
            in.get("static_condenstation", static_condenstation);
            in.get("stol", stol_);

            // SQP
            Scalar_t upper_bound_value = 1;
            in.get("upper_bound", upper_bound_value);

            in.require("from", [&](Input &node) {
                node.require("space", *space_from);
                auto fun_from = std::make_shared<FEModelFunction_t>(space_from);
                node.require("function", *fun_from);
                functions_.push_back(fun_from);

                node.get("solver", *linear_solver_);
            });

            in.require("to", [&](Input &node) {
                node.require("space", *space_to);
                auto fun_to = std::make_shared<FEModelFunction_t>(space_to);
                node.require("function", *fun_to);
                functions_.push_back(fun_to);

                node.get("solver", *qp_solver_);
            });

            in.get("transfer", transfer_);

            if (!transfer_.init(space_from, space_to)) {
                Utopia::Abort("FETransfer::init failed!");
            }

            // SQP
            auto upper_bound = std::make_shared<Vector_t>();
            functions_[1]->space()->create_vector(*upper_bound);
            upper_bound->set(upper_bound_value);
            box_.upper_bound() = upper_bound;
            qp_solver_->set_box_constraints(box_);

            active_set_.init(layout(*upper_bound));
            active_set_.verbose(true);
        }

        bool run() {
            auto &&comm = functions_[0]->space()->comm();

            std::size_t n = functions_.size();
            std::vector<SubSystem> systems(n);

            for (std::size_t i = 0; i < n; ++i) {
                auto &s = systems[i];
                auto &f = functions_[i];

                f->create_solution_vector(s.solution);
                f->hessian_and_gradient(s.solution, s.hessian, s.g);
                f->space()->apply_zero_constraints(s.g);
            }

            Vector_t transferred(layout(systems[1].g), 0.);
            Vector_t residual(layout(systems[0].g), 0.);
            Vector_t correction(layout(systems[0].g), 0.);

            Vector_t to_residual;
            Vector_t to_correction;
            Vector_t solution_diff;

            ////////////////////////////////////////////////////
            auto &to_system = systems[1];
            auto &from_system = systems[0];

            Scalar_t prev_norm = 1e6;

            IO<FunctionSpace> io(*functions_[0]->space());
            io.set_output_path("iterations.e");

            ////////////////////////////////////////////////////
            /// 1) Prepare first condensed linear system matrix
            //////////////////////////////////////////////////
            // Jacobian/Hessian of coupled problem
            Matrix_t jacobian;

            // remove identity from rows
            functions_[1]->space()->apply_constraints(to_system.hessian, 0);

            transfer_.apply(to_system.hessian, jacobian);

            jacobian += from_system.hessian;

            linear_solver_->update(make_ref(jacobian));

            // reset identity to rows
            functions_[1]->space()->apply_constraints(to_system.hessian, 1);
            ////////////////////////////////////////////////////

            // Intialize QPSolver
            qp_solver_->update(make_ref(to_system.hessian));

            for (int k = 0; k < picard_iterations; ++k) {
                comm.root_print("---------------------------------------\n");
                //////////////////////////////////////////////////
                // 2) Solve condensed system (for correction)
                //////////////////////////////////////////////////

                functions_[0]->gradient(from_system.solution, from_system.g);

                // contribution from system 1
                transfer_.apply_transpose(to_system.g, residual);
                residual += from_system.g;

                // Negative gradient!
                residual = -residual;

                {
                    Scalar_t norm_residual = norm2(residual);
                    comm.root_print(std::to_string(k) + ") norm_residual =" + std::to_string(norm_residual));
                }

                linear_solver_->apply(residual, correction);

                // Update solution with correction
                from_system.solution += correction;

                Scalar_t norm_correction = norm2(correction);
                comm.root_print(std::to_string(k) + ") norm_correction =" + std::to_string(norm_correction));

                if (norm_correction < stol_) {
                    comm.root_print("Converged!");
                    break;
                }

                // Write to disk
                io.write(from_system.solution, k, k);

                //////////////////////////////////////////////////
                // 3) Transfer to qp problem
                //////////////////////////////////////////////////
                transfer_.apply(from_system.solution, transferred);

                //////////////////////////////////////////////////
                // 4) Solving for x (not correction)
                //////////////////////////////////////////////////

                to_residual = to_system.hessian * transferred;
                to_system.solution = transferred;

                functions_[1]->space()->apply_constraints(to_residual);
                qp_solver_->apply(to_residual, to_system.solution);

                solution_diff = transferred - to_system.solution;

                {
                    Scalar_t norm_solution_diff = norm2(solution_diff);
                    comm.root_print(std::to_string(k) + ") norm_solution_diff =" + std::to_string(norm_solution_diff));
                }

                //////////////////////////////////////////////////
                // 5) Compute gradient of coupled system
                //////////////////////////////////////////////////

                // Constrained stress
                functions_[1]->gradient(to_system.solution, to_system.g);
                // functions_[1]->gradient(transferred, to_residual);

                // to_system.g = to_system.hessian * solution_diff;
                // to_system.g = to_residual - to_system.g;

                {
                    Scalar_t norm_g_1 = norm2(to_system.g);
                    comm.root_print(std::to_string(k) + ") norm_g_1 =" + std::to_string(norm_g_1));
                }
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
        std::shared_ptr<OmniLinearSolver_t> linear_solver_;
        std::shared_ptr<OmniQPSolver_t> qp_solver_;

        // SQP
        BoxConstraints<Vector_t> box_;
        ActiveSet<Vector_t> active_set_;

        int picard_iterations{2};
        bool first_order{false};
        Scalar_t dumping{1};
        bool write_matlab{false};
        bool augment_diagonal{true};
        bool static_condenstation{false};
        Scalar_t stol_{1e-6};
    };
}  // namespace utopia

#endif  // UTOPIA_SADDLE_POINT_QP_SOLVE_HPP
