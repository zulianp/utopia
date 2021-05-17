#ifndef UTOPIA_TIME_DEPENDENT_NONLINEAR_PROBLEM_HPP
#define UTOPIA_TIME_DEPENDENT_NONLINEAR_PROBLEM_HPP

#include "utopia_Input.hpp"
#include "utopia_ui.hpp"

#include "utopia_Field.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_fe_Environment.hpp"

#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_MatrixTransformer.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_ProblemBase.hpp"

namespace utopia {

    // values have to be organized to fit a Newton solver
    // f(x) + b; -> J(x) * s = -(f(x) + b), x = x + s

    template <class FunctionSpace>
    class TimeDependentNonlinearProblem : public TimeDependentProblem<FunctionSpace> {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using OmniAssembler_t = utopia::OmniAssembler<FunctionSpace>;
        using LinearSolver_t = utopia::LinearSolver<Matrix_t, Vector_t>;
        using OmniLinearSolver_t = utopia::OmniLinearSolver<Matrix_t, Vector_t>;
        using Environment_t = utopia::Environment<FunctionSpace>;
        using IO_t = utopia::IO<FunctionSpace>;
        using Super = utopia::TimeDependentProblem<FunctionSpace>;

        using ImplicitEulerIntegrator_t = utopia::ImplicitEulerIntegrator<FunctionSpace>;
        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;

        void read(Input &in) override {
            Super::read(in);

            std::string integration_scheme = "ImplicitEuler";
            node.get("integration_scheme", integration_scheme);

            assert(this->space());

            if (integration_scheme == "Newmark") {
                integrator_ = std::make_shared<NewmarkIntegrator_t>(this->space());
            } else
            // if (integration_scheme == "ImplicitEuler")
            {
                integrator_ = std::make_shared<ImplicitEulerIntegrator_t>(this->space());
            }

            output_path_ = this->output_dir() + ("/" + this->name() + ".e");
            io_->set_output_path(output_path_);
            in.get("solver", *linear_solver_);
            in.get("verbose", verbose_);
        }

        void set_environment(const std::shared_ptr<Environment_t> &env) override { integrator_->set_environment(env); }

        bool init() override {
            if (!Super::init()) return false;
            assert(this->delta_time() > 0);
            return true;
        }

        bool assemble_material() {
            assert(Scalar_t(norm2(*this->solution())) == 0.0);
            if (!integrator_->assemble(*this->solution(), *this->jacobian(), *this->fun())) {
                return false;
            }

            *this->fun() *= -1.0;
            return true;
        }

        bool assemble_operators() override { return true; }

        inline bool must_update_system() const { return this->is_first_time_step() || system_modified_flag_; }

        bool prepare_system() override {
            compute_system();
            compute_residual();
            return true;
        }

        bool apply_constraints() override {
            if (must_update_system()) {
                assert(!empty(system_));
                this->space()->apply_constraints(system_);
            }

            assert(!empty(residual_));
            this->space()->apply_zero_constraints(residual_);
            return true;
        }

        bool update() override { return assemble_operators() && prepare_system() && apply_constraints() && solve(); }

        void apply_transformers() {
            for (auto &trafo : transformers) {
                trafo->apply(*this->jacobian());
            }
        }

        void condensed_system_built() override { system_modified_flag_ = true; }

        bool solve() override {
            increment_.set(0.0);

            if (must_update_system()) {
                linear_solver_->update(make_ref(system_));
                system_modified_flag_ = false;
            }

            bool ok = linear_solver_->apply(residual_, increment_);
            assert(ok);

            (*this->solution()) += increment_;

            if (verbose_) {
                utopia::out() << "Step:\t" << this->current_time().step() << "\tTime:\t" << this->current_time().get()
                              << "\tSum increment " << Scalar_t(sum(increment_)) << '\n';
            }

            return ok;
        }

        bool export_result() const override {
            if (this->integrate_all_before_output()) {
                if (this->complete()) {
                    return io_->write(*this->solution());
                }
            } else {
                return io_->write(*this->solution(), this->current_time().step(), this->current_time().get());
            }

            return true;
        }

        inline std::vector<std::shared_ptr<Matrix_t>> operators() override {
            std::vector<std::shared_ptr<Matrix_t>> ret{this->jacobian(), mass_matrix()};
            return ret;
        }

        inline bool is_linear() const override { return true; }
        inline bool has_transformers() const { return !transformers.empty(); }

        TimeDependentNonlinearProblem(const std::shared_ptr<FunctionSpace> &space)
            : Super(space),
              linear_solver_(std::make_shared<OmniLinearSolver_t>()),
              io_(std::make_shared<IO_t>(*space)) {}

    private:
        std::shared_ptr<TimeDependentFunction_t> integrator_;
        std::shared_ptr<LinearSolver_t> linear_solver_;
        std::shared_ptr<IO_t> io_;
        std::vector<std::unique_ptr<MatrixTransformer<Matrix_t>>> transformers;

        Path output_path_;
        bool verbose_{false};
        bool system_modified_flag_{false};
    };

}  // namespace utopia
#endif  // UTOPIA_TIME_DEPENDENT_NONLINEAR_PROBLEM_HPP