#ifndef UTOPIA_BARRIER_MULTIGRID_HPP
#define UTOPIA_BARRIER_MULTIGRID_HPP
#include "utopia_ConvergenceReason.hpp"
#include "utopia_Core.hpp"

#include "utopia_IPTransfer.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Level.hpp"
#include "utopia_LinearMultiLevel.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_Recorder.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_Utils.hpp"

#include "utopia_LineSearchBoxProjection.hpp"
#include "utopia_LogBarrierFunctionBase.hpp"
#include "utopia_LogBarrierFunctionFactory.hpp"

#include <cassert>
#include <ctime>

namespace utopia {
    /**
     * @brief      The class for  BarrierMultigrid solver.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class BarrierMultigrid final : public NonLinearSolver<Vector>,
                                   public VariableBoundSolverInterface<Vector>,
                                   public Clonable {
        using Layout = typename Traits<Vector>::Layout;
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        using Super = utopia::QPSolver<Matrix, Vector>;

        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;
        using Smoother = utopia::IterativeSolver<Matrix, Vector>;
        using SmootherPtr = std::shared_ptr<Smoother>;
        using IterativeSolver = utopia::IterativeSolver<Matrix, Vector>;
        using Level = utopia::Level<Matrix, Vector>;
        using Transfer = utopia::Transfer<Matrix, Vector>;
        using IPTransfer = utopia::IPTransfer<Matrix, Vector>;
        using VariableBoundSolverInterface = utopia::VariableBoundSolverInterface<Vector>;

        using LogBarrierBase = utopia::LogBarrierBase<Matrix, Vector>;

        // class NonLinearSmoother : public Configurable {
        //     void update(const std::shared_ptr<const Matrix> &op) {}
        //     bool apply(const Vector &rhs, Vector &x) { return false; }

        //     void read(Input &in) override {}
        // };

        // using DefaultSmoother = NonLinearSmoother;

        class LevelMemory {
        public:
            Matrix matrix;
            Vector residual;
            Vector correction;
        };

    public:
        BarrierMultigrid(
            const std::shared_ptr<LinearSolver> &coarse_solver
            // , const std::shared_ptr<NonLinearSmoother> &nl_smoother = std::make_shared<DefaultSmoother>()
            )
            : coarse_solver_(coarse_solver)
        // , nl_smoother_(nl_smoother)
        {}

        ~BarrierMultigrid() override = default;

        void read(Input &in) override {
            // Super::read(in);

            std::string barrier_function_type;

            Options()
                .add_option("barrier_function_type",
                            barrier_function_type,
                            "Type of LogBarrier. Options={LogBarrierFunctionWithSelection|LogBarrierFunction}")
                .add_option("debug", debug_, "TODO")
                .add_option("use_coarse_space", use_coarse_space_, "TODO")
                .add_option("pre_smoothing_steps", pre_smoothing_steps_, "TODO")
                .add_option("post_smoothing_steps", post_smoothing_steps_, "TODO")
                .add_option("use_non_linear_residual", use_non_linear_residual_, "TODO")
                .parse(in);

            barrier_ = LogBarrierFunctionFactory<Matrix, Vector>::new_log_barrier(barrier_function_type);
            barrier_->read(in);

            // if (nl_smoother_) {
            //     in.get("nl_smoother", *nl_smoother_);
            // }
        }

        void set_transfer_operators(const std::vector<std::shared_ptr<Transfer>> &transfer_operators) {
            this->transfer_operators_ = transfer_operators;
        }

        // void init(const std::shared_ptr<NonLinearSmoother> &nl_smoother) { nl_smoother_ = nl_smoother; }

        // void update(const std::shared_ptr<const Matrix> &op) override { Super::update(op); }

        bool solve(Function<Matrix, Vector> &fun, Vector &x) {
            bool converged = false;

            // std::string mg_header_message = "BarrierMultigrid: " + std::to_string(n_levels()) + " levels";
            // this->init_solver(mg_header_message, {" it. ", "|| P(x-g)-x ||", "|| u - u_old ||", "# active sets"});

            // Vector x_old = x;
            // Vector grad;
            // Scalar g_norm = this->criticality_measure_inf(x, grad);

            // if (this->verbose()) {
            //     PrintInfo::print_iter_status(0, {g_norm, 0});
            // }

            ensure_defaults();

            auto &&box = this->get_box_constraints();

            barrier_->set_box_constraints(make_ref(box));
            barrier_->reset();

            line_search_projection_->set_box_constraints(make_ref(box));

            const int n_coarse_levels = this->n_levels() - 1;

            Vector correction;

            //  Function quantities
            Matrix f_hessian;
            Vector f_gradient;
            Vector f_diag_hessian;

            // Barrier quantities
            Vector b_gradient;
            Matrix b_hessian;
            b_hessian.identity(square_matrix_layout(layout(x)), 0.);

            Vector b_diag_hessian;

            Vector diag_op;
            Vector residual;
            Vector delta_x;

            Matrix f_hessian_coarse;
            Vector b_hessian_coarse;

            for (SizeType it = 1; it < this->max_it() && !converged; ++it) {
                fun.hessian(x, f_hessian);

                if (debug_) {
                    write("X_" + std::to_string(it) + ".m", x);
                }

                transfer_operators_[n_coarse_levels - 1]->restrict(f_hessian, f_hessian_coarse);

                if (!use_non_linear_residual_) {
                    fun.gradient(x, f_gradient);
                }

                f_diag_hessian = diag(f_hessian);

                for (int s = 0; s < mg_steps_; ++s) {
                    ////////////////////////////////////////////////////////////
                    // Pre-smoothing
                    ////////////////////////////////////////////////////////////
                    for (int ps = 0; ps < pre_smoothing_steps(); ++ps) {
                        if (use_non_linear_residual_) {
                            fun.gradient(x, f_gradient);
                        }

                        residual = f_gradient;
                        barrier_->gradient(x, residual);
                        residual *= -1;

                        // FIXME
                        barrier_->hessian(x, b_hessian);
                        b_diag_hessian = diag(b_hessian);

                        diag_op = f_diag_hessian + b_diag_hessian;
                        diag_op = 1. / diag_op;

                        correction = e_mul(diag_op, residual);

                        Scalar alpha = line_search_projection_->compute(x, correction);

                        x += alpha * correction;

                        if (!use_non_linear_residual_) {
                            f_gradient -= alpha * (f_hessian * correction);
                        }
                    }

                    ////////////////////////////////////////////////////////////
                    // Coarse grid correction
                    ////////////////////////////////////////////////////////////
                    if (use_coarse_space_) {
                        b_hessian *= 0.;
                        barrier_->hessian(x, b_hessian);

                        // FIXME
                        b_diag_hessian = diag(b_hessian);

                        auto &&mem = memory_[n_coarse_levels - 1];
                        auto &&transfer = transfer_operators_[n_coarse_levels - 1];

                        if (use_non_linear_residual_) {
                            fun.gradient(x, f_gradient);
                        }

                        residual = -f_gradient;
                        barrier_->gradient(x, residual);

                        transfer->restrict(residual, mem.residual);
                        transfer->restrict(b_diag_hessian, b_hessian_coarse);

                        mem.matrix = f_hessian_coarse;
                        mem.matrix += diag(b_hessian_coarse);
                        mem.correction.zeros(layout(mem.residual));

                        coarse_solver_->solve(mem.matrix, mem.residual, mem.correction);

                        transfer->interpolate(mem.correction, correction);

                        // FIXME

                        delta_x = correction;

                        if (box.has_upper_bound()) {
                            delta_x = x + correction;
                            delta_x = min(delta_x, *box.upper_bound());
                            delta_x -= x;
                        }

                        if (box.has_lower_bound()) {
                            delta_x += x;
                            delta_x = max(delta_x, *box.lower_bound());
                            delta_x -= x;
                        }

                        Scalar alpha = line_search_projection_->compute(x, correction);

                        delta_x *= dumping_ * line_search_projection_weight_;
                        delta_x += ((1 - line_search_projection_weight_) * alpha) * correction;

                        x += delta_x;

                        if (!use_non_linear_residual_) {
                            f_gradient -= (f_hessian * delta_x);
                        }
                    }

                    ////////////////////////////////////////////////////////////
                    // Post-smoothing
                    ////////////////////////////////////////////////////////////
                    // for (int ps = 0; ps < post_smoothing_steps(); ++ps) {
                    //     if (use_non_linear_residual_) {
                    //         fun.gradient(x, f_gradient);
                    //     }

                    //     residual = -f_gradient;
                    //     barrier_->gradient(x, residual);
                    // }
                }

                // converged = this->check_convergence(it, g_norm, 1, 1);

                // if (!converged) {
                //     x_old = x;
                // }
            }

            return converged;
        }

        void init_memory(const Layout &) override { ensure_memory(); }

        // inline bool step(Function<Matrix, Vector> &fun, Vector &x) {
        //     // nl_smoother_->set_box_constraints(this->get_box_constraints());
        //     // nl_smoother_->sweeps(pre_smoothing_steps());
        //     // nl_smoother_->smooth(rhs, x);

        //     // TODO

        //     // nl_smoother_->sweeps(post_smoothing_steps());
        //     // nl_smoother_->smooth(rhs, x);
        //     return true;
        // }

        void pre_smoothing_steps(const SizeType n) { pre_smoothing_steps(n); }
        void post_smoothing_steps(const SizeType n) { post_smoothing_steps(n); }

        int n_levels() const { return transfer_operators_.size() + 1; }

    public:
        BarrierMultigrid *clone() const override {
            auto ret = utopia::make_unique<BarrierMultigrid>(std::shared_ptr<LinearSolver>(coarse_solver_->clone()));
            return ret.release();
        }

        /**
         * @brief      Function changes direct solver needed for coarse grid solve.
         *
         * @param[in]  linear_solver  The linear solver.
         *
         * @return
         */
        bool change_coarse_solver(const std::shared_ptr<LinearSolver> &coarse_solver) {
            return coarse_solver_ = (coarse_solver);
        }

        Scalar criticality_measure_inf(const Vector &x, const Vector &g) {
            Vector help = x - g;

            this->project(*this->constraints_.lower_bound(), *this->constraints_.upper_bound(), help);
            help -= x;

            return norm2(help);
        }

        inline int pre_smoothing_steps() const { return pre_smoothing_steps_; }
        inline int post_smoothing_steps() const { return post_smoothing_steps_; }

    protected:
        std::shared_ptr<LinearSolver> coarse_solver_;
        // std::shared_ptr<NonLinearSmoother> nl_smoother_;
        std::vector<std::shared_ptr<Transfer>> transfer_operators_;
        std::vector<LevelMemory> memory_;
        std::shared_ptr<LogBarrierBase> barrier_;
        std::shared_ptr<LineSearchBoxProjection<Vector>> line_search_projection_;

        bool debug_{false};
        bool use_coarse_space_{true};
        bool use_non_linear_residual_{false};

        int mg_steps_{1};
        int pre_smoothing_steps_{3};
        int post_smoothing_steps_{3};

        Scalar dumping_{0.98};
        Scalar line_search_projection_weight_{0.9};

        void ensure_defaults() {
            if (!barrier_) {
                barrier_ = std::make_shared<LogBarrier<Matrix, Vector>>();
            }

            if (!line_search_projection_) {
                line_search_projection_ = std::make_shared<LineSearchBoxProjection<Vector>>();
            }

            ensure_memory();
        }

        void ensure_memory() {
            if (memory_.size() != transfer_operators_.size()) {
                memory_.resize(transfer_operators_.size());
            }
        }

    };  // namespace utopia

}  // namespace utopia

#endif  // UTOPIA_BARRIER_MULTIGRID_HPP
