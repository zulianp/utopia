#ifndef UTOPIA_TWO_FIELD_SPIN_HPP
#define UTOPIA_TWO_FIELD_SPIN_HPP

#include <iomanip>
#include <limits>
#include <utopia.hpp>
#include "utopia_Backtracking.hpp"
#include "utopia_Newton.hpp"

namespace utopia {

    class SolutionStatusSPIN {
    public:
        int nl_iterates_field1{0};
        int nl_iterates_field2{0};

        int num_linear_solves_field1{0};
        int num_linear_solves_field2{0};

        int sum_linear_its_field1{0};
        int sum_linear_its_field2{0};

        double cond_number_precond{0.0};
        double cond_number_unprecond{0.0};

        SolutionStatusSPIN(SolutionStatus &sol_status) : global_(sol_status) {}

        void clear() { global_.clear(); }

        inline void describe(std::ostream &os) const {
            global_.describe(os);

            os << "nl_iterates_field1:       " << nl_iterates_field1 << "\n";
            os << "nl_iterates_field2:       " << nl_iterates_field2 << "\n";
            os << "num_linear_solves_field1:        " << num_linear_solves_field1 << "  \n";
            os << "num_linear_solves_field2: " << num_linear_solves_field2 << "\n";
            os << "sum_linear_its_field1:  " << sum_linear_its_field1 << "\n";
            os << "sum_linear_its_field2: " << sum_linear_its_field2 << " \n";
        }

        SolutionStatus &global_;
    };

    template <class Matrix, class Vector = typename Traits<Matrix>::Vector>
    class TwoFieldSPIN final : public NonLinearSolver<Vector>, public InexactNewtonInterface<Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        using LinearSolver = utopia::MatrixFreeLinearSolver<Vector>;
        using NonLinearSolver = utopia::NonLinearSolver<Vector>;

        using FieldNonLinearSolver = utopia::Newton<Matrix, Vector>;

        using LSStrategy = utopia::Backtracking<Vector>;

        typedef utopia::ExtendedFunction<Matrix, Vector> Fun;
        using FunPtr = std::shared_ptr<Fun>;

        using Communicator = typename Traits<Vector>::Communicator;

    public:
        // to be moved around
        class FunctionOperator final : public Operator<Vector> {
        public:
            FunctionOperator(TwoFieldSPIN &parent, const std::function<void(const Vector &, Vector &)> operator_action)
                : parent(parent), operator_action_(operator_action) {}

            bool apply(const Vector &rhs, Vector &ret) const override {
                operator_action_(rhs, ret);
                return true;
            }

            inline Communicator &comm() override { return parent.comm(); }
            inline const Communicator &comm() const override { return parent.comm(); }

            inline Size size() const override { return parent.size(); }
            inline Size local_size() const override { return parent.local_size(); }

        private:
            TwoFieldSPIN &parent;
            std::function<void(const Vector &, Vector &)> operator_action_;
        };

        class FunctionPreconditionerSPIN final : public Preconditioner<Vector> {
        public:
            FunctionPreconditionerSPIN(TwoFieldSPIN &parent,
                                       const std::function<void(const Vector &, Vector &)> operator_action)
                : parent(parent), operator_action_(operator_action) {}

            bool apply(const Vector &rhs, Vector &ret) override {
                operator_action_(rhs, ret);
                return true;
            }

            FunctionPreconditionerSPIN *clone() const override { return new FunctionPreconditionerSPIN(*this); }

        private:
            TwoFieldSPIN &parent;
            std::function<void(const Vector &, Vector &)> operator_action_;
        };

        TwoFieldSPIN(const std::shared_ptr<LinearSolver> &linear_solver =
                         std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>(),
                     const std::shared_ptr<FieldNonLinearSolver> &nonlinear_solver_field1 =
                         std::make_shared<Newton<Matrix, Vector>>(),
                     const std::shared_ptr<FieldNonLinearSolver> &nonlinear_solver_field2 =
                         std::make_shared<Newton<Matrix, Vector>>())
            : NonLinearSolver(),
              InexactNewtonInterface<Vector>(),
              alpha_(1.0),
              additive_(true),
              global_linear_solver_(std::move(linear_solver)),
              nonlinear_solver_field1_(std::move(nonlinear_solver_field1)),
              nonlinear_solver_field2_(std::move(nonlinear_solver_field2)),
              sol_status_spin_(this->solution_status_global()) {}

        bool solve(Function<Matrix, Vector> &fun_global, Vector &x_global) {
            using namespace utopia;

            x_field1_ = fun_field1_->initial_guess();
            x_field2_ = fun_field2_->initial_guess();

            init_memory(x_global, x_field1_, x_field2_);
            Scalar g_norm, s_norm, r_norm, g0_norm, objective_val, objective_val_new_;
            Scalar norm_lin_resid, g_norm_old, nom, denom, c_max;
            SizeType it = 0;

            bool SPIN_flg = false;
            bool converged = false;

            fun_global.update(x_global);
            fun_global.value_and_gradient(x_global, objective_val, grad_global_);

            norms2(grad_global_, step_global_, g_norm, s_norm);
            g0_norm = g_norm;
            r_norm = 1.0;

            auto SPIN_operator = build_SPIN_operator();

            this->init_solver(
                "SPIN",
                {" it. ", "|| g ||", "   J    ", "r_norm", "|| s_global || ", "alpha_k", "eta", "rho", "SPIN_flg"});

            if (this->verbose_) {
                PrintInfo::print_iter_status(it, {g_norm, objective_val, r_norm, s_norm, alpha_});
            }

            it++;

            if (use_switch_) {
                linear_solver_global2_ = std::make_shared<GMRES<Matrix, Vector>>();

                if (this->use_block_precond_global()) {
                    auto precond = this->build_linear_preconditioner_operator();
                    linear_solver_global2_->set_preconditioner(precond);
                } else {
                    linear_solver_global2_->pc_type("hypre");
                }

                linear_solver_global2_->verbose(false);
                linear_solver_global2_->atol(1e-12);
                linear_solver_global2_->max_it(1000);

                c_max = norm_infty(x_field2_);

                if (c_max < c_max_tol_) {
                    eta_ = 0.0;
                    rho_ = 1.0;
                }
            }

            while (!converged) {
                if (!use_switch_ or (!((eta_ < eps1_) || (rho_ > eps2_)))) {
                    SPIN_flg = true;
                    global_to_field1_(x_global, x_field1_);
                    global_to_field2_(x_global, x_field2_);

                    if (additive_) {
                        s1_ = x_field1_;  // initial guess for nonlinear_solve_field1
                        nonlinear_solve_field1();
                        s2_ = x_field2_;  // initial guess for nonlinear_solve_field2

                        nonlinear_solve_field2();
                        field1_to_global_(s1_, g_spin_global_);
                        field2_to_global_(s2_, g_spin_global_);
                        // g_spin_global_ = [u*; c*]

                    } else  // multiplicative
                    {
                        s1_ = x_field1_;
                        nonlinear_solve_field1();

                        g_spin_global_ = x_global;
                        field1_to_global_(s1_, g_spin_global_);
                        global_to_field2_(g_spin_global_, s2_);

                        // s2_ = x_field2_;
                        nonlinear_solve_field2();
                        field2_to_global_(s2_, g_spin_global_);
                    }

                    // evaluating Jacobians/hesssians
                    if (this->use_exact_hessian()) {
                        fun_global.hessian(g_spin_global_, H_global_);
                        fun_field1_->hessian(s1_, H_field1_);
                        fun_field2_->hessian(s2_, H_field2_);
                    } else {
                        fun_global.hessian(x_global, H_global_);
                        fun_field1_->hessian(x_field1_, H_field1_);
                        fun_field2_->hessian(x_field2_, H_field2_);
                    }

                    // evaluate nonlinearly preconditioned residual
                    g_spin_global_ -= x_global;

                    // zeroing search direction / initial guess
                    step_global_.set(0.0);

                    // set verbosity
                    if (this->verbosity_level() > VERBOSITY_LEVEL_NORMAL && this->verbose() == true) {
                        if (auto *ls_solver =
                                dynamic_cast<IterativeSolver<Matrix, Vector> *>(global_linear_solver_.get())) {
                            ls_solver->verbose(true);
                        }
                    }

                    // inexact newton stuff
                    if (this->has_forcing_strategy()) {
                        if (auto *iterative_solver =
                                dynamic_cast<IterativeSolver<Matrix, Vector> *>(global_linear_solver_.get())) {
                            global_ls_atol_estimate_ = this->estimate_ls_atol(g_norm, it);
                            iterative_solver->atol(global_ls_atol_estimate_);
                            // iterative_solver->verbose(true);
                        }
                    }

                    // std::cout << "SPIN solve .... \n";
                    global_linear_solver_->solve(*SPIN_operator, g_spin_global_, step_global_);

                    this->solution_status_.num_linear_solves++;
                    if (auto *it_solver =
                            dynamic_cast<IterativeSolver<Matrix, Vector> *>(global_linear_solver_.get())) {
                        auto sol_status_ls = it_solver->solution_status();
                        // std::cout << "sol_status_ls.iterates  " << sol_status_ls.iterates
                        //           << "  \n";
                        this->solution_status_.sum_linear_its += sol_status_ls.iterates;
                    }

                    if (ls_strategy_) {
                        // bool ls_flg = ls_strategy_->get_alpha(fun_global, grad_global_,
                        //                                       x_global, step_global_,
                        //                                       alpha_);

                        bool ls_flg = ls_strategy_->get_alpha(fun_global,
                                                              grad_global_,
                                                              x_global,
                                                              step_global_,
                                                              objective_val,
                                                              objective_val_new_,
                                                              g_help_,
                                                              alpha_);

                        if (ls_flg == false) {
                            std::cout << "taking AM step \n";
                            step_global_ = g_spin_global_;

                            bool ls_flg = ls_strategy_->get_alpha(fun_global,
                                                                  grad_global_,
                                                                  x_global,
                                                                  step_global_,
                                                                  objective_val,
                                                                  objective_val_new_,
                                                                  g_help_,
                                                                  alpha_);
                        }

                        if (ls_flg == false) {
                            if (!linear_solver_global2_) {
                                linear_solver_global2_ = std::make_shared<GMRES<Matrix, Vector>>();
                                // linear_solver_global2_->pc_type("hypre");
                                linear_solver_global2_->verbose(false);
                                linear_solver_global2_->atol(1e-12);
                                linear_solver_global2_->max_it(100000);
                            }
                            std::cout << "Taking Newton step \n";

                            // TODO:: recompute in case of ASPEN
                            // if (this->use_block_precond_global()) {
                            //   // in order to setup block preconditioner, ideally could be
                            //   // extracted from other stuff
                            //   global_to_field1_(x_global, x_field1_);
                            //   global_to_field2_(x_global, x_field2_);
                            //   fun_field1_->hessian(x_field1_,
                            //   H_field1_);
                            //   fun_field2_->hessian(x_field2_,
                            //   H_field2_);
                            // }

                            g_spin_global_ = -1.0 * grad_global_;

                            // linear_solver_global2_->verbose(false);
                            // std::cout << "Newton solve .... \n";
                            linear_solver_global2_->solve(H_global_, g_spin_global_, step_global_);

                            // ls_flg = ls_strategy_->get_alpha(fun_global, grad_global_,
                            // x_global,
                            //                                  step_global_, alpha_);

                            bool ls_flg = ls_strategy_->get_alpha(fun_global,
                                                                  grad_global_,
                                                                  x_global,
                                                                  step_global_,
                                                                  objective_val,
                                                                  objective_val_new_,
                                                                  g_help_,
                                                                  alpha_);
                        }

                        if (ls_flg == false) {
                            std::cout << "Taking gradient step \n";
                            step_global_ = -1.0 * grad_global_;

                            // bool ls_flg = ls_strategy_->get_alpha(
                            //     fun_global, grad_global_, x_global, step_global_, alpha_);

                            bool ls_flg = ls_strategy_->get_alpha(fun_global,
                                                                  grad_global_,
                                                                  x_global,
                                                                  step_global_,
                                                                  objective_val,
                                                                  objective_val_new_,
                                                                  g_help_,
                                                                  alpha_);
                        }

                        step_global_ *= alpha_;
                        x_global += step_global_;
                    } else {
                        // update
                        if (fabs(alpha_ - 1.0) < std::numeric_limits<Scalar>::epsilon()) {
                            x_global += step_global_;
                        } else {
                            step_global_ *= alpha_;
                            x_global += step_global_;
                        }
                    }
                } else {
                    SPIN_flg = false;
                    fun_global.hessian(x_global, H_global_);

                    if (this->use_block_precond_global()) {
                        // in order to setup block preconditioner, ideally could be extracted
                        // from other stuff
                        global_to_field1_(x_global, x_field1_);
                        global_to_field2_(x_global, x_field2_);
                        fun_field1_->hessian(x_field1_, H_field1_);
                        fun_field2_->hessian(x_field2_, H_field2_);
                    }

                    step_global_.set(0.0);

                    // inexact newton stuff
                    if (this->has_forcing_strategy()) {
                        if (auto *iterative_solver =
                                dynamic_cast<IterativeSolver<Matrix, Vector> *>(linear_solver_global2_.get())) {
                            global_ls_atol_estimate_ = this->estimate_ls_atol(g_norm, it);
                            iterative_solver->atol(global_ls_atol_estimate_);
                        }
                    }

                    // to be changed
                    g_spin_global_ = -1.0 * grad_global_;

                    linear_solver_global2_->solve(H_global_, g_spin_global_, step_global_);

                    this->solution_status_.num_linear_solves++;
                    if (auto *it_solver =
                            dynamic_cast<IterativeSolver<Matrix, Vector> *>(linear_solver_global2_.get())) {
                        auto sol_status_ls = it_solver->solution_status();
                        this->solution_status_.sum_linear_its += sol_status_ls.iterates;
                    }

                    if (ls_strategy_) {
                        // ls_strategy_->get_alpha(fun_global, grad_global_, x_global,
                        //                         step_global_, alpha_);

                        bool ls_flg = ls_strategy_->get_alpha(fun_global,
                                                              grad_global_,
                                                              x_global,
                                                              step_global_,
                                                              objective_val,
                                                              objective_val_new_,
                                                              g_help_,
                                                              alpha_);

                        step_global_ *= alpha_;
                        x_global += step_global_;
                    } else {
                        // update x
                        if (fabs(alpha_ - 1.0) < std::numeric_limits<Scalar>::epsilon()) {
                            x_global += step_global_;
                        } else {
                            step_global_ *= alpha_;
                            x_global += step_global_;
                        }
                    }
                }  // end of inexact Newton's stuff

                if (use_switch_) {
                    // reusing g_spin_global_ vector to save memory
                    g_spin_global_ = H_global_ * (1. / alpha_ * step_global_);
                    g_spin_global_ = grad_global_ + g_spin_global_;
                    g_norm_old = g_norm;
                }

                fun_global.update(x_global);

                if (ls_strategy_) {
                    objective_val = objective_val_new_;
                    grad_global_ = g_help_;
                } else {
                    fun_global.value_and_gradient(x_global, objective_val, grad_global_);
                }

                if (use_switch_) {
                    norms2(grad_global_, step_global_, g_spin_global_, g_norm, s_norm, norm_lin_resid);
                } else {
                    norms2(grad_global_, step_global_, g_norm, s_norm);
                }
                r_norm = g_norm / g0_norm;

                if (use_switch_) {
                    global_to_field2_(x_global, x_field2_);
                    c_max = norm_infty(x_field2_);

                    if (c_max < c_max_tol_) {
                        eta_ = 0.0;
                        rho_ = 1.0;
                    } else {
                        nom = std::abs(g_norm - norm_lin_resid);
                        denom = g_norm_old;

                        eta_ = nom / denom;
                        rho_ = (g_norm_old - g_norm) / (g_norm_old - norm_lin_resid);
                    }
                }

                if (this->verbose_) {
                    PrintInfo::print_iter_status(
                        it, {g_norm, objective_val, r_norm, s_norm, alpha_, eta_, rho_, static_cast<Scalar>(SPIN_flg)});
                }

                converged = this->check_convergence(it, g_norm, r_norm, s_norm);
                if (alpha_ < alpha_tol_) {
                    converged = true;
                    this->solution_status_.iterates = it;
                    this->solution_status_.gradient_norm = g_norm;
                    this->solution_status_.relative_gradient_norm = r_norm;
                    this->solution_status_.step_norm = s_norm;
                    this->exit_solver(it, 7);
                }

                it++;
            }

            return true;
        }

        void read(Input &in) override {
            NonLinearSolver::read(in);
            in.get("dumping", alpha_);
            in.get("additive", additive_);
            in.get("use_switch", use_switch_);

            if (ls_strategy_) {
                in.get("line-search", *ls_strategy_);
            }

            if (global_linear_solver_) {
                in.get("global_linear_solver", *global_linear_solver_);
            }

            if (nonlinear_solver_field1_) {
                in.get("nonlinear_solver_field1", *nonlinear_solver_field1_);
            }

            if (nonlinear_solver_field2_) {
                in.get("nonlinear_solver_field2", *nonlinear_solver_field2_);
            }
        }

        void print_usage(std::ostream &os) const override {
            NonLinearSolver::print_usage(os);

            this->print_param_usage(os, "dumping", "real", "Default step size.", "1.0");
            this->print_param_usage(os, "aditive", "bool", "Type of preconditioner (additive/multiplicative)", "true");

            this->print_param_usage(os,
                                    "use_switch",
                                    "bool",
                                    "Use adaptive switch between SPIN precond. and "
                                    "pure inexact Newton's method",
                                    "false");

            this->print_param_usage(os, "line-search", "LSStrategy", "Input parameters for line-search strategy.", "-");
            this->print_param_usage(os,
                                    "global_linear_solver",
                                    "MatrixFreeLinearSolver",
                                    "Input parameters for global linear solver.",
                                    "-");
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
        }

        void set_line_search_strategy(const std::shared_ptr<LSStrategy> &strategy) { ls_strategy_ = strategy; }

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

        virtual Size size() const { return size_; }

        virtual Size local_size() const { return local_size_; }

    private:
        void init_memory(const Vector &x_global, const Vector &x_field1, const Vector &x_field2) {
            const Layout layout_global = layout(x_global);
            const Layout layout_field1 = layout(x_field1);
            const Layout layout_field2 = layout(x_field2);

            // init of linear solver
            // NonLinearSolver::init_memory(layout_global);

            // init global quantities
            step_global_.zeros(layout_global);
            grad_global_.zeros(layout_global);
            g_spin_global_.zeros(layout_global);
            g_help_.zeros(layout_global);

            if (global_linear_solver_) global_linear_solver_->init_memory(layout_global);

            if (ls_strategy_) ls_strategy_->init_memory(layout_global);

            x_field1_.zeros(layout_field1);
            s1_.zeros(layout_field1);
            if (nonlinear_solver_field1_) nonlinear_solver_field1_->init_memory(layout_field1);

            x_field2_.zeros(layout_field2);
            frac_marker_local_.zeros(layout_field2);
            s2_.zeros(layout_field2);

            if (nonlinear_solver_field2_) nonlinear_solver_field2_->init_memory(layout_field2);

            comm_ = std::shared_ptr<Communicator>(x_global.comm().clone());

            size_.set_dims(2);
            size_.set(0, x_global.size());
            size_.set(1, x_global.size());

            local_size_.set_dims(2);
            local_size_.set(0, x_global.local_size());
            local_size_.set(1, x_global.local_size());
        }

        std::shared_ptr<FunctionOperator> build_SPIN_operator() {
            if (additive_) {
                return build_additive_SPIN_operator();
            } else {
                return build_multiplicative_SPIN_operator();
            }
        }

        std::shared_ptr<FunctionPreconditionerSPIN> build_linear_preconditioner_operator() {
            if (additive_) {
                return build_additive_linear_precond();
            } else {
                return build_multiplicative_linear_precond();
            }
        }

        std::shared_ptr<FunctionOperator> build_additive_SPIN_operator() {
            std::function<void(const Vector &, Vector &)> my_func_additive = [this](const Vector &p, Vector &result) {
                result = H_global_ * p;

                this->global_to_field1_(result, s1_);
                this->global_to_field2_(result, s2_);

                x_field1_.set(0.0);
                x_field2_.set(0.0);

                // additive version
                auto ls_field1 = nonlinear_solver_field1_->linear_solver();
                auto ls_field2 = nonlinear_solver_field2_->linear_solver();

                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field1.get())) {
                    ls_solver->atol(1e-1 * global_ls_atol_estimate_);
                    ls_solver->rtol(action_rtol_);
                    if (this->verbosity_level() == VERBOSITY_LEVEL_DEBUG && this->verbose() == true) {
                        ls_solver->verbose(true);
                        out() << "--- field 1 - linear solve in SPIN mat_apply ---- \n";
                    }
                }

                ls_field1->solve(H_field1_, s1_, x_field1_);

                sol_status_spin_.num_linear_solves_field1 += 1;
                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field1.get())) {
                    const auto &sol_status_ls = ls_solver->solution_status();
                    sol_status_spin_.sum_linear_its_field1 += sol_status_ls.iterates;
                }

                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field2.get())) {
                    if (this->verbosity_level() == VERBOSITY_LEVEL_DEBUG && this->verbose() == true) {
                        ls_solver->verbose(true);
                        out() << "--- field 2 - linear solve in SPIN mat_apply ---- \n";
                    }

                    ls_solver->atol(1e-1 * global_ls_atol_estimate_);
                    ls_solver->rtol(action_rtol_);
                }

                ls_field2->solve(H_field2_, s2_, x_field2_);

                sol_status_spin_.num_linear_solves_field2 += 1;
                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field2.get())) {
                    const auto &sol_status_ls = ls_solver->solution_status();
                    sol_status_spin_.sum_linear_its_field2 += sol_status_ls.iterates;
                }

                this->field1_to_global_(x_field1_, result);
                this->field2_to_global_(x_field2_, result);
            };

            return std::make_shared<FunctionOperator>(*this, my_func_additive);
        }

        std::shared_ptr<FunctionOperator> build_multiplicative_SPIN_operator() {
            std::function<void(const Vector &, Vector &)> my_func_additive = [this](const Vector &p, Vector &result) {
                result = H_global_ * p;

                this->global_to_field1_(result, s1_);
                this->global_to_field2_(result, s2_);

                x_field1_.set(0.0);

                // additive version
                auto ls_field1 = nonlinear_solver_field1_->linear_solver();
                auto ls_field2 = nonlinear_solver_field2_->linear_solver();

                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field1.get())) {
                    ls_solver->atol(1e-1 * global_ls_atol_estimate_);
                    ls_solver->rtol(action_rtol_);
                    if (this->verbosity_level() == VERBOSITY_LEVEL_DEBUG && this->verbose() == true) {
                        ls_solver->verbose(true);
                        out() << "--- field 1 - linear solve in SPIN mat_apply ---- \n";
                    }
                }

                ls_field1->solve(H_field1_, s1_, x_field1_);

                sol_status_spin_.num_linear_solves_field1 += 1;
                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field1.get())) {
                    const auto &sol_status_ls = ls_solver->solution_status();
                    sol_status_spin_.sum_linear_its_field1 += sol_status_ls.iterates;
                }

                result.set(0.0);
                this->field1_to_global_(x_field1_, result);
                result = H_global_ * result;
                this->global_to_field2_(result, x_field2_);
                s2_ -= x_field2_;

                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field2.get())) {
                    if (this->verbosity_level() == VERBOSITY_LEVEL_DEBUG && this->verbose() == true) {
                    }
                    ls_solver->atol(1e-1 * global_ls_atol_estimate_);
                    ls_solver->rtol(action_rtol_);
                }

                x_field2_.set(0.0);

                ls_field2->solve(H_field2_, s2_, x_field2_);

                sol_status_spin_.num_linear_solves_field2 += 1;
                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field2.get())) {
                    const auto &sol_status_ls = ls_solver->solution_status();
                    sol_status_spin_.sum_linear_its_field2 += sol_status_ls.iterates;
                }

                this->field1_to_global_(x_field1_, result);
                this->field2_to_global_(x_field2_, result);
            };

            return std::make_shared<FunctionOperator>(*this, my_func_additive);
        }

        std::shared_ptr<FunctionPreconditionerSPIN> build_additive_linear_precond() {
            std::function<void(const Vector &, Vector &)> my_func_additive_precond = [this](const Vector &rhs,
                                                                                            Vector &result) {
                this->global_to_field1_(rhs, s1_);
                this->global_to_field2_(rhs, s2_);

                x_field1_.set(0.0);
                x_field2_.set(0.0);

                // additive version
                auto ls_field1 = nonlinear_solver_field1_->linear_solver();
                auto ls_field2 = nonlinear_solver_field2_->linear_solver();

                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field1.get())) {
                    ls_solver->atol(1e-1 * global_ls_atol_estimate_);
                    ls_solver->rtol(action_rtol_);
                    if (this->verbosity_level() == VERBOSITY_LEVEL_DEBUG && this->verbose() == true) {
                        ls_solver->verbose(true);
                        out() << "--- field 1 - linear solve in SPIN mat_apply ---- \n";
                    }
                }

                ls_field1->solve(H_field1_, s1_, x_field1_);

                sol_status_spin_.num_linear_solves_field1 += 1;
                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field1.get())) {
                    const auto &sol_status_ls = ls_solver->solution_status();
                    sol_status_spin_.sum_linear_its_field1 += sol_status_ls.iterates;
                }

                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field2.get())) {
                    if (this->verbosity_level() == VERBOSITY_LEVEL_DEBUG && this->verbose() == true) {
                        ls_solver->verbose(true);
                        out() << "--- field 2 - linear solve in SPIN mat_apply ---- \n";
                    }

                    ls_solver->atol(1e-1 * global_ls_atol_estimate_);
                    ls_solver->rtol(action_rtol_);
                }

                ls_field2->solve(H_field2_, s2_, x_field2_);

                sol_status_spin_.num_linear_solves_field2 += 1;
                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field2.get())) {
                    const auto &sol_status_ls = ls_solver->solution_status();
                    sol_status_spin_.sum_linear_its_field2 += sol_status_ls.iterates;
                }

                this->field1_to_global_(x_field1_, result);
                this->field2_to_global_(x_field2_, result);
            };

            return std::make_shared<FunctionPreconditionerSPIN>(*this, my_func_additive_precond);
        }

        std::shared_ptr<FunctionPreconditionerSPIN> build_multiplicative_linear_precond() {
            std::function<void(const Vector &, Vector &)> my_func_mult = [this](const Vector &rhs, Vector &result) {
                // result = H_global_ * p;

                this->global_to_field1_(rhs, s1_);
                this->global_to_field2_(rhs, s2_);

                x_field1_.set(0.0);

                // additive version
                auto ls_field1 = nonlinear_solver_field1_->linear_solver();
                auto ls_field2 = nonlinear_solver_field2_->linear_solver();

                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field1.get())) {
                    ls_solver->atol(1e-1 * global_ls_atol_estimate_);
                    ls_solver->rtol(action_rtol_);
                    if (this->verbosity_level() == VERBOSITY_LEVEL_DEBUG && this->verbose() == true) {
                        ls_solver->verbose(true);
                        out() << "--- field 1 - linear solve in SPIN mat_apply ---- \n";
                    }
                }

                ls_field1->solve(H_field1_, s1_, x_field1_);

                sol_status_spin_.num_linear_solves_field1 += 1;
                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field1.get())) {
                    const auto &sol_status_ls = ls_solver->solution_status();
                    sol_status_spin_.sum_linear_its_field1 += sol_status_ls.iterates;
                }

                result.set(0.0);
                this->field1_to_global_(x_field1_, result);
                result = H_global_ * result;
                this->global_to_field2_(result, x_field2_);
                s2_ -= x_field2_;

                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field2.get())) {
                    if (this->verbosity_level() == VERBOSITY_LEVEL_DEBUG && this->verbose() == true) {
                        ls_solver->verbose(true);
                        out() << "--- field 2 - linear solve in SPIN mat_apply ---- \n";
                    }
                    ls_solver->atol(1e-1 * global_ls_atol_estimate_);
                    ls_solver->rtol(action_rtol_);
                }

                x_field2_.set(0.0);

                ls_field2->solve(H_field2_, s2_, x_field2_);

                sol_status_spin_.num_linear_solves_field2 += 1;
                if (auto *ls_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(ls_field2.get())) {
                    const auto &sol_status_ls = ls_solver->solution_status();
                    sol_status_spin_.sum_linear_its_field2 += sol_status_ls.iterates;
                }

                this->field1_to_global_(x_field1_, result);
                this->field2_to_global_(x_field2_, result);
            };

            return std::make_shared<FunctionPreconditionerSPIN>(*this, my_func_mult);
        }

        void nonlinear_solve_field1() {
            if (this->verbosity_level() > VERBOSITY_LEVEL_NORMAL && this->verbose() == true) {
                out() << "------ nonlinear solve for the first field ----- \n";
                nonlinear_solver_field1_->verbose(true);
            }
            nonlinear_solver_field1_->solve(*fun_field1_, s1_);

            const auto &sol_status = nonlinear_solver_field1_->solution_status();
            // sol_status.describe(std::cout);

            sol_status_spin_.nl_iterates_field1 += sol_status.iterates;
            sol_status_spin_.num_linear_solves_field1 += sol_status.num_linear_solves;
            sol_status_spin_.sum_linear_its_field1 += sol_status.sum_linear_its;
        }

        void nonlinear_solve_field2() {
            if (this->verbosity_level() > VERBOSITY_LEVEL_NORMAL && this->verbose() == true) {
                out() << "------ nonlinear solve for the second field ----- \n ";
                nonlinear_solver_field2_->verbose(true);
            }
            nonlinear_solver_field2_->solve(*fun_field2_, s2_);

            const auto &sol_status = nonlinear_solver_field2_->solution_status();
            // sol_status.describe(std::cout);

            sol_status_spin_.nl_iterates_field2 += sol_status.iterates;
            sol_status_spin_.num_linear_solves_field2 += sol_status.num_linear_solves;
            sol_status_spin_.sum_linear_its_field2 += sol_status.sum_linear_its;
        }

    public:
        virtual Communicator &comm() {
            assert(comm_);
            return *comm_;
        }

        virtual const Communicator &comm() const {
            assert(comm_);
            return *comm_;
        }

        virtual VerbosityLevel verbosity_level() const { return _verbosity_level; }

        virtual void verbosity_level(const VerbosityLevel &level) { _verbosity_level = level; }

        void action_rtol(const Scalar &action_rtol) { action_rtol_ = action_rtol; }
        Scalar action_rtol() const { return action_rtol_; }

        void alpha_tol(const Scalar &alpha_tol) { alpha_tol_ = alpha_tol; }
        Scalar alpha_tol() const { return alpha_tol_; }

        void additive_precond(const bool &additive) { additive_ = additive; }
        bool additive_precond() const { return additive_; }

        void multiplicative_precond(const bool &multiplicative) { additive_ = !multiplicative; }
        bool multiplicative_precond() const { return !additive_; }

        void use_switch(const bool &use_switch) { use_switch_ = use_switch; }
        bool use_switch() const { return use_switch_; }

        void use_exact_hessian(const bool &use_exact_hessian) { use_exact_hessian_ = use_exact_hessian; }
        bool use_exact_hessian() const { return use_exact_hessian_; }

        void use_block_precond_global(const bool &use_block_precond_global) {
            use_block_precond_global_ = use_block_precond_global;
        }
        bool use_block_precond_global() const { return use_block_precond_global_; }

        SolutionStatus &solution_status_global() { return this->solution_status_; }
        const SolutionStatusSPIN &solution_status() { return sol_status_spin_; }

    private:
        Scalar alpha_;
        bool additive_;
        bool use_switch_{false};
        bool use_exact_hessian_{false};
        bool use_block_precond_global_{false};

        Scalar action_rtol_{1e-14};
        Scalar alpha_tol_{1e-6};
        Scalar global_ls_atol_estimate_{1e-14};
        Scalar c_max_tol_{0.98};

        Scalar eta_{1.};
        Scalar rho_{0.0};
        Scalar eps1_{0.2};
        Scalar eps2_{0.8};

        std::shared_ptr<LinearSolver> global_linear_solver_;
        std::shared_ptr<GMRES<Matrix, Vector>> linear_solver_global2_;

        std::shared_ptr<FieldNonLinearSolver> nonlinear_solver_field1_;
        std::shared_ptr<FieldNonLinearSolver> nonlinear_solver_field2_;

        std::shared_ptr<LSStrategy> ls_strategy_;

        FunPtr fun_field1_;
        FunPtr fun_field2_;

        std::function<void(const Vector &, Vector &)> field1_to_global_;
        std::function<void(const Vector &, Vector &)> field2_to_global_;
        std::function<void(const Vector &, Vector &)> global_to_field1_;
        std::function<void(const Vector &, Vector &)> global_to_field2_;

        Vector step_global_, grad_global_, g_spin_global_, g_help_;

        Vector x_field1_, s1_;
        Vector x_field2_, s2_, frac_marker_local_;

        Matrix H_global_, H_field1_, H_field2_;
        Matrix D_inv_field1_, D_inv_field2_;

        Size size_, local_size_;
        std::shared_ptr<Communicator> comm_;

        VerbosityLevel _verbosity_level{VERBOSITY_LEVEL_NORMAL};

        SolutionStatusSPIN sol_status_spin_;
    };

}  // namespace utopia

#endif  // UTOPIA_TWO_FIELD_SPIN_HPP
