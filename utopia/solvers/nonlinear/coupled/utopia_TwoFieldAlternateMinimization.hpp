#ifndef UTOPIA_ALTERNATE_MINIMIZATION_HPP
#define UTOPIA_ALTERNATE_MINIMIZATION_HPP

#include <iomanip>
#include <limits>
#include <utopia.hpp>
#include "utopia_TwoFieldSPIN.hpp"

#include "utopia_Traits.hpp"

namespace utopia {

    template <class Matrix, class Vector = typename Traits<Matrix>::Vector>
    class TwoFieldAlternateMinimization final : public NonLinearSolver<Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        using LinearSolver = utopia::MatrixFreeLinearSolver<Vector>;
        using NonLinearSolver = utopia::NonLinearSolver<Vector>;

        using FieldNonLinearSolver = utopia::NewtonBase<Matrix, Vector>;

        using LSStrategy = utopia::LSStrategy<Vector>;

        typedef utopia::Function<Matrix, Vector> Fun;
        using FunPtr = std::shared_ptr<Fun>;

        using NonLinearSolver::check_convergence;

    public:
        TwoFieldAlternateMinimization(const std::shared_ptr<FieldNonLinearSolver> &nonlinear_solver_field1 =
                                          std::make_shared<Newton<Matrix, Vector>>(),
                                      const std::shared_ptr<FieldNonLinearSolver> &nonlinear_solver_field2 =
                                          std::make_shared<Newton<Matrix, Vector>>())
            : NonLinearSolver(),
              nonlinear_solver_field1_(std::move(nonlinear_solver_field1)),
              nonlinear_solver_field2_(std::move(nonlinear_solver_field2)),
              sol_status_spin_(this->solution_status_global()) {}

        bool solve(Function<Matrix, Vector> &fun_global, Vector &x_global) {
            using namespace utopia;

            // x_field1_ = fun_field1_->initial_guess();
            // x_field2_ = fun_field2_->initial_guess();

            fun_field1_->create_vector(x_field1_);
            fun_field2_->create_vector(x_field2_);

            init_memory(x_global, x_field1_, x_field2_);
            Scalar g_norm, r_norm, g0_norm, objective_val, diff_field_1, diff_field_2;

            SizeType it = 0;
            bool converged = false;

            fun_global.update(x_global);
            fun_global.gradient(x_global, grad_global_);
            fun_global.value(x_global, objective_val);

            g_norm = norm2(grad_global_);
            g0_norm = g_norm;
            r_norm = 1.0;

            this->init_solver("AlternateMinimization",
                              {" it. ", "|| g ||", "   J    ", "r_norm", "||disp_diff ||_{inf} ", "||c ||_{inf}"});

            if (this->verbose_) {
                PrintInfo::print_iter_status(it, {g_norm, objective_val, r_norm});
            }

            it++;

            while (!converged) {
                global_to_field1_(x_global, x_field1_);
                x_field1_old_ = x_field1_;
                nonlinear_solve_field1();

                field1_to_global_(x_field1_, x_global);
                global_to_field2_(x_global, x_field2_);
                x_field2_old_ = x_field2_;

                nonlinear_solve_field2();
                field2_to_global_(x_field2_, x_global);

                fun_global.update(x_global);
                fun_global.value(x_global, objective_val);
                fun_global.gradient(x_global, grad_global_);

                g_norm = norm2(grad_global_);
                r_norm = g_norm / g0_norm;

                x_field1_old_ = x_field1_ - x_field1_old_;
                x_field2_old_ = x_field2_ - x_field2_old_;

                diff_field_1 = norm_infty(x_field1_old_);
                diff_field_2 = norm_infty(x_field2_old_);

                if (this->verbose_) {
                    PrintInfo::print_iter_status(it, {g_norm, objective_val, r_norm, diff_field_1, diff_field_2});
                }

                converged = this->check_convergence(it, g_norm, r_norm, 1.0, diff_field_1, diff_field_2);

                it++;
            }

            return true;
        }

        void read(Input &in) override {
            NonLinearSolver::read(in);

            in.get("field1_diff_tol", field1_diff_tol_);
            in.get("field2_diff_tol", field2_diff_tol_);

            if (nonlinear_solver_field1_) {
                in.get("nonlinear_solver_field1", *nonlinear_solver_field1_);
            }

            if (nonlinear_solver_field2_) {
                in.get("nonlinear_solver_field2", *nonlinear_solver_field2_);
            }
        }

        void print_usage(std::ostream &os) const override {
            NonLinearSolver::print_usage(os);

            this->print_param_usage(os,
                                    "nonlinear_solver_field1",
                                    "NonLinearSolver",
                                    "Input parameters for nonlinear solver used for first field",
                                    "-");
            this->print_param_usage(os,
                                    "nonlinear_solver_field1",
                                    "NonLinearSolver",
                                    "Input parameters for nonlinear solver used for first field",
                                    "-");

            this->print_param_usage(os,
                                    "field1_diff_tol",
                                    "Scalar",
                                    "Tolerance used for measuring difference in abs. value of fields",
                                    "1e-14");

            this->print_param_usage(os,
                                    "field2_diff_tol",
                                    "Scalar",
                                    "Tolerance used for measuring difference in abs. value of fields",
                                    "1e-14");
        }

        void set_field_functions(const FunPtr &field_1, const FunPtr &field_2) {
            fun_field1_ = field_1;
            fun_field2_ = field_2;
        }

        void set_transfers(const std::function<void(const Vector &, Vector &)> &field1_to_global,
                           const std::function<void(const Vector &, Vector &)> &field2_to_global,
                           const std::function<void(const Vector &, Vector &)> &global_to_field1,
                           const std::function<void(const Vector &, Vector &)> &global_to_field2) {
            field1_to_global_ = field1_to_global;
            field2_to_global_ = field2_to_global;
            global_to_field1_ = global_to_field1;
            global_to_field2_ = global_to_field2;
        }

        virtual VerbosityLevel verbosity_level() const { return _verbosity_level; }

        virtual void verbosity_level(const VerbosityLevel &level) { _verbosity_level = level; }

    private:
        void init_memory(const Vector &x_global, const Vector &x_field1, const Vector &x_field2) {
            const Layout layout_global = layout(x_global);
            const Layout layout_field1 = layout(x_field1);
            const Layout layout_field2 = layout(x_field2);

            // global
            grad_global_.zeros(layout_global);

            x_field1_.zeros(layout_field1);
            x_field1_old_.zeros(layout_field1);
            if (nonlinear_solver_field1_) nonlinear_solver_field1_->init_memory(layout_field1);

            x_field2_.zeros(layout_field2);
            x_field2_old_.zeros(layout_field2);

            if (nonlinear_solver_field2_) nonlinear_solver_field2_->init_memory(layout_field2);
        }

        void nonlinear_solve_field1() {
            if (this->verbosity_level() > VERBOSITY_LEVEL_NORMAL && this->verbose() == true && mpi_world_rank() == 0) {
                out() << "------ nonlinear solve for the first field ----- \n";
                nonlinear_solver_field1_->verbose(true);
            }
            nonlinear_solver_field1_->solve(*fun_field1_, x_field1_);

            const auto &sol_status = nonlinear_solver_field1_->solution_status();
            // sol_status.describe(std::cout);

            sol_status_spin_.nl_iterates_field1 += sol_status.iterates;
            sol_status_spin_.num_linear_solves_field1 += sol_status.num_linear_solves;
            sol_status_spin_.sum_linear_its_field1 += sol_status.sum_linear_its;
        }

        void nonlinear_solve_field2() {
            if (this->verbosity_level() > VERBOSITY_LEVEL_NORMAL && this->verbose() == true && mpi_world_rank() == 0) {
                out() << "------ nonlinear solve for the second field ----- \n ";
                nonlinear_solver_field2_->verbose(true);
            }
            nonlinear_solver_field2_->solve(*fun_field2_, x_field2_);

            const auto &sol_status = nonlinear_solver_field2_->solution_status();
            // sol_status.describe(std::cout);

            sol_status_spin_.nl_iterates_field2 += sol_status.iterates;
            sol_status_spin_.num_linear_solves_field2 += sol_status.num_linear_solves;
            sol_status_spin_.sum_linear_its_field2 += sol_status.sum_linear_its;
        }

        bool check_convergence(const SizeType &it,
                               const Scalar &g_norm,
                               const Scalar &r_norm,
                               const Scalar &s_norm,
                               const Scalar &diff_field_1,
                               const Scalar &diff_field_2) {
            bool converged = NonLinearSolver::check_convergence(it, g_norm, r_norm, 1.0);

            if (!converged and ((diff_field_1 <= field1_diff_tol_) or (diff_field_2 <= field2_diff_tol_))) {
                converged = true;
                NonLinearSolver::check_convergence(it, g_norm, r_norm, 1e-15);
            }

            return converged;
        }

    public:
        SolutionStatus &solution_status_global() { return this->solution_status_; }
        const SolutionStatusSPIN &solution_status() { return sol_status_spin_; }

        void field1_diff_tol(const Scalar &field1_diff_tol) { field1_diff_tol_ = field1_diff_tol; }
        void field2_diff_tol(const Scalar &field2_diff_tol) { field2_diff_tol_ = field2_diff_tol; }

        Scalar field1_diff_tol() const { return field1_diff_tol_; }
        Scalar field2_diff_tol() const { return field2_diff_tol_; }

    private:
        std::shared_ptr<FieldNonLinearSolver> nonlinear_solver_field1_;
        std::shared_ptr<FieldNonLinearSolver> nonlinear_solver_field2_;

        FunPtr fun_field1_;
        FunPtr fun_field2_;

        std::function<void(const Vector &, Vector &)> field1_to_global_;
        std::function<void(const Vector &, Vector &)> field2_to_global_;
        std::function<void(const Vector &, Vector &)> global_to_field1_;
        std::function<void(const Vector &, Vector &)> global_to_field2_;

        Vector grad_global_;
        Vector x_field1_, x_field1_old_, x_field2_, x_field2_old_;

        SolutionStatusSPIN sol_status_spin_;
        VerbosityLevel _verbosity_level{VERBOSITY_LEVEL_NORMAL};

        Scalar field1_diff_tol_{1e-20};
        Scalar field2_diff_tol_{1e-20};
    };

}  // namespace utopia

#endif
