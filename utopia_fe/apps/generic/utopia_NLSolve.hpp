#ifndef UTOPIA_NL_SOLVE_HPP
#define UTOPIA_NL_SOLVE_HPP

// Base and Algebra
#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia.hpp"

// UtopiaFE
#include "utopia_FEModelFunction.hpp"
#include "utopia_fe_Environment.hpp"

// #include "utopia_OmniLinearSolver.hpp"

namespace utopia {

    template <class FunctionSpace>
    class NLSolve : public Configurable {
    public:
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using FEFunctionInterface_t = utopia::FEFunctionInterface<FunctionSpace>;

        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;

        using NewtonBase_t = utopia::NewtonBase<Matrix_t, Vector_t>;
        using Newton_t = utopia::Newton<Matrix_t, Vector_t>;

        inline std::shared_ptr<Environment_t> &environment() {
            assert(environment_);
            return environment_;
        }
        inline const Communicator_t &comm() const { return comm_; }
        inline Communicator_t &comm() { return comm_; }

        void init(const std::shared_ptr<FEFunctionInterface_t> &function) { function_ = function; }

        void read(Input &in) override {
            in.get("solver", *solver_);
            in.get("verbose", verbose_);
            in.get("use_pseudo_newton", use_pseudo_newton_);
        }

        inline void set_solver(const std::shared_ptr<NewtonBase_t> &solver) { solver_ = solver; }

        bool solve_trivial() {
            Vector_t x;
            function_->create_solution_vector(x);
            // zero out boundary conditions
            x.set(0.0);

            this->status("Assemblying linear problem");

            Matrix_t A;
            Vector_t rhs;
            bool ok = function_->hessian_and_gradient(x, A, rhs);
            assert(ok);

            rhs *= -1;
            // function_->space()->write("rhs.e", rhs);

            // apply boundary conditions to rhs
            function_->space()->apply_constraints(rhs);

            this->status("Solving linear problem");
            // Solve linear problem
            ok = solver_->linear_solver()->solve(A, rhs, x);
            assert(ok);

            this->status("Reporting solution");

            function_->report_solution(x);
            return ok;
        }

        bool solve_trivial_pseudo_newton() {
            Vector_t x;
            function_->create_solution_vector(x);

            this->status("Assemblying linear problem");

            Matrix_t H;  // Hessian/Jacobian
            Vector_t g;  // Postive gradient / Residual
            bool ok = function_->hessian_and_gradient(x, H, g);
            assert(ok);

            function_->space()->write("rhs.e", g);

            Vector_t c(layout(x), 0.0);  // Correction

            this->status("Solving linear problem (pseudo-Newton)");
            // Solve linear problem
            ok = solver_->linear_solver()->solve(H, g, c);
            assert(ok);

            // Correct the solution with respect to the negative gradient
            x -= c;

            this->status("Reporting solution");

            function_->report_solution(x);
            return ok;
        }

        bool solve() {
            if (function_->is_linear() && !function_->is_time_dependent()) {
                // Tivial problem, lets keep it simple
                if (use_pseudo_newton_) {
                    return solve_trivial_pseudo_newton();
                } else {
                    return solve_trivial();
                }
            } else {
                Vector_t x;
                function_->create_solution_vector(x);

                function_->setup_IVP(x);

                this->status("Solving nonlinear problem");

                int n_time_steps = 0;

                do {
                    this->status("Timestep: " + std::to_string(n_time_steps++));

                    if (!solver_->solve(*function_, x)) {
                        error("Solver failed to solve");
                        return false;
                    }

                    function_->update_IVP(x);
                    function_->report_solution(x);
                } while (!function_->is_IVP_solved());
            }

            return true;
        }

        NLSolve() : environment_(std::make_shared<Environment_t>()) { init_defaults(); }

        inline bool verbose() const { return verbose_; }
        inline void verbose(const bool val) { verbose_ = val; }
        virtual std::string name() const { return "NLSolve"; }

        void status(const std::string &message) const {
            if (verbose_) {
                this->comm().root_print("[Status] " + name() + ": " + message, utopia::out().stream());
            }
        }

        void warning(const std::string &message) const {
            if (verbose_) {
                this->comm().root_print("[Warning] " + name() + ": " + message, utopia::err().stream());
            }
        }

        void error(const std::string &message) const {
            if (verbose_) {
                this->comm().root_print("[Error] " + name() + ": " + message, utopia::err().stream());
            }
        }

    private:
        Communicator_t comm_;
        std::shared_ptr<Environment_t> environment_;
        std::shared_ptr<FEFunctionInterface_t> function_;
        std::shared_ptr<NewtonBase_t> solver_;
        bool verbose_{true};
        bool use_pseudo_newton_{false};

        void init_defaults() {
            auto linear_solver = std::make_shared<OmniLinearSolver_t>();
            auto newton = std::make_shared<Newton_t>(linear_solver);
            solver_ = newton;
        }
    };
}  // namespace utopia

#endif