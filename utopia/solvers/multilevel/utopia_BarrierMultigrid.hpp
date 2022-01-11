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

        using Super = utopia::NonLinearSolver<Vector>;

        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;
        using Smoother = utopia::IterativeSolver<Matrix, Vector>;
        using SmootherPtr = std::shared_ptr<Smoother>;
        using IterativeSolver = utopia::IterativeSolver<Matrix, Vector>;
        using Level = utopia::Level<Matrix, Vector>;
        using Transfer = utopia::Transfer<Matrix, Vector>;
        using IPTransfer = utopia::IPTransfer<Matrix, Vector>;
        using VariableBoundSolverInterface = utopia::VariableBoundSolverInterface<Vector>;

        using LogBarrierBase = utopia::LogBarrierBase<Matrix, Vector>;

        class LevelMemory {
        public:
            Matrix matrix;
            Vector residual;
            Vector correction;
        };

    public:
        BarrierMultigrid(const std::shared_ptr<LinearSolver> &coarse_solver) : coarse_solver_(coarse_solver) {}

        ~BarrierMultigrid() override = default;

        void read(Input &in) override {
            Super::read(in);

            std::string barrier_function_type = "BoundedLogBarrier";
            // std::string barrier_function_type = "LogBarrier";

            Options()
                .add_option("barrier_function_type",
                            barrier_function_type,
                            "Type of LogBarrier. Options={LogBarrierWithSelection|LogBarrier}")
                .add_option("debug", debug_, "Enable/Disable debug ouput.")
                .add_option("use_coarse_space", use_coarse_space_, "Use the coarse space for multigrid acceleration.")
                .add_option("pre_smoothing_steps", pre_smoothing_steps_, "Number of pre-smoothing  steps.")
                .add_option("post_smoothing_steps", post_smoothing_steps_, "Number of post-smoothing  steps.")
                .add_option("mg_steps", mg_steps_, "Number of multigrid iterations per Newton step.")
                .add_option("use_non_linear_residual",
                            use_non_linear_residual_,
                            "Use non-linear residual, instead of linear within the multigrid iteration.")
                .parse(in);

            barrier_ = LogBarrierFunctionFactory<Matrix, Vector>::new_log_barrier(barrier_function_type);
            barrier_->read(in);
        }

        void set_transfer_operators(const std::vector<std::shared_ptr<Transfer>> &transfer_operators) {
            this->transfer_operators_ = transfer_operators;
        }
        bool solve(Function<Matrix, Vector> &fun, Vector &x) {
            bool converged = false;

            std::string mg_header_message = "BarrierMultigrid: " + std::to_string(n_levels()) + " levels";
            this->init_solver(mg_header_message, {" it. ", "|| g ||"});

            ensure_defaults();

            auto &&box = this->get_box_constraints();

            barrier_->set_box_constraints(make_ref(box));
            barrier_->reset();

            line_search_projection_->set_box_constraints(make_ref(box));

            const int n_coarse_levels = this->n_levels() - 1;

            IterationState state;

            //  Function quantities
            auto &f_gradient = state.f_gradient;
            auto &f_hessian = state.f_hessian;
            auto &f_diag_hessian = state.f_diag_hessian;

            // Barrier quantities
            auto &b_gradient = state.b_gradient;

            auto &residual = state.residual;
            auto &f_hessian_coarse = state.f_hessian_coarse;

            for (SizeType it = 1; it < this->max_it() && !converged; ++it) {
                fun.hessian(x, f_hessian);
                transfer_operators_[n_coarse_levels - 1]->restrict(f_hessian, f_hessian_coarse);

                if (!use_non_linear_residual_) {
                    fun.gradient(x, f_gradient);
                }

                f_diag_hessian = diag(f_hessian);

                for (int s = 0; s < mg_steps_; ++s) {
                    ////////////////////////////////////////////////////////////
                    // Pre-smoothing
                    ////////////////////////////////////////////////////////////
                    this->smooth(fun, x, state, pre_smoothing_steps());

                    ////////////////////////////////////////////////////////////
                    // Coarse grid correction
                    ////////////////////////////////////////////////////////////
                    if (use_coarse_space_) {
                        coarse_correction(fun, x, state);
                    }

                    ////////////////////////////////////////////////////////////
                    // Post-smoothing
                    ////////////////////////////////////////////////////////////
                    this->smooth(fun, x, state, post_smoothing_steps());
                }

                barrier_->update_barrier();

                fun.gradient(x, f_gradient);

                b_gradient.set(0.);
                barrier_->gradient(x, b_gradient);

                residual = f_gradient + b_gradient;

                Scalar g_norm = norm2(residual);
                PrintInfo::print_iter_status(it, {g_norm});

                converged = this->check_convergence(it, g_norm, 1, 1);
            }

            if (this->verbose()) {
                std::stringstream ss;
                state.alpha_stats.describe(ss);

                x.comm().root_print(ss.str());
            }

            if (debug_) {
                rename("r", b_gradient);

                write("R.m", b_gradient);

                rename("g", f_gradient);
                f_gradient *= -1;
                write("G.m", f_gradient);
            }

            return converged;
        }

        void init_memory(const Layout &) override { ensure_memory(); }

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

        UTOPIA_NVCC_PRIVATE

        std::shared_ptr<LinearSolver> coarse_solver_;
        std::vector<std::shared_ptr<Transfer>> transfer_operators_;
        std::vector<LevelMemory> memory_;
        std::shared_ptr<LogBarrierBase> barrier_;
        std::shared_ptr<LineSearchBoxProjection<Vector>> line_search_projection_;

        bool debug_{true};
        bool use_coarse_space_{false};
        bool use_non_linear_residual_{false};

        int mg_steps_{1};
        int pre_smoothing_steps_{3};
        int post_smoothing_steps_{3};

        Scalar dumping_{0.98};
        Scalar line_search_projection_weight_{0.9};

        class AlphaStats {
        public:
            Scalar min_{1};
            Scalar max_{0};
            Scalar avg_{0};

            SizeType n_searches_{0};

            void reset() {
                min_ = 1;
                max_ = 0;
                avg_ = 0;
                n_searches_ = 0;
            }

            void alpha_found(const Scalar alpha) {
                Scalar scale = n_searches_ / (n_searches_ + 1.);
                min_ = std::min(alpha, min_);
                max_ = std::max(alpha, max_);
                avg_ = scale * avg_ + alpha / (n_searches_ + 1);
                n_searches_++;
            }

            void describe(std::ostream &os) const {
                os << "Alpha) min: " << min_ << ", max: " << max_ << ", avg: " << avg_ << ", n: " << n_searches_
                   << "\n";
            }
        };

        void ensure_defaults() {
            if (!barrier_) {
                barrier_ = std::make_shared<BoundedLogBarrier<Matrix, Vector>>();
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

        class IterationState {
        public:
            Vector correction;

            //  Function quantities
            Matrix f_hessian;
            Vector f_gradient;
            Vector f_diag_hessian;

            // Barrier quantities
            Vector b_gradient;
            Vector b_diag_hessian;

            Vector diag_op;
            Vector residual;
            Vector delta_x;

            Matrix f_hessian_coarse;
            Vector b_hessian_coarse;

            AlphaStats alpha_stats;
        };

        void handle_boundary(Vector &v) const {
            // FIXME
            Write<Vector> w(v);
            auto r = range(v);
            if (r.inside(0)) {
                v.set(0, 0);
            }

            if (r.inside(v.size() - 1)) {
                v.set(v.size() - 1, 0);
            }
        }

        void coarse_correction(Function<Matrix, Vector> &fun, Vector &x, IterationState &state) {
            auto &f_gradient = state.f_gradient;
            auto &f_hessian = state.f_hessian;

            // auto &b_gradient = state.b_gradient;
            auto &b_diag_hessian = state.b_diag_hessian;

            auto &correction = state.correction;
            auto &residual = state.residual;

            auto &delta_x = state.delta_x;
            auto &f_hessian_coarse = state.f_hessian_coarse;
            auto &b_hessian_coarse = state.b_hessian_coarse;

            auto &&box = this->get_box_constraints();

            const int n_coarse_levels = this->n_levels() - 1;

            b_diag_hessian *= 0.;
            barrier_->hessian_diag(x, b_diag_hessian);

            // FIXME
            handle_boundary(b_diag_hessian);

            auto &&mem = memory_[n_coarse_levels - 1];
            auto &&transfer = transfer_operators_[n_coarse_levels - 1];

            if (use_non_linear_residual_) {
                fun.gradient(x, f_gradient);
            }

            residual = f_gradient;
            barrier_->gradient(x, residual);
            residual *= -1;

            handle_boundary(residual);

            transfer->restrict(residual, mem.residual);
            transfer->restrict(b_diag_hessian, b_hessian_coarse);

            mem.matrix = f_hessian_coarse;
            mem.matrix += diag(b_hessian_coarse);
            mem.correction.zeros(layout(mem.residual));

            coarse_solver_->solve(mem.matrix, mem.residual, mem.correction);

            transfer->interpolate(mem.correction, correction);

            handle_boundary(correction);

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
                f_gradient += (f_hessian * delta_x);
            }
        }

        void smooth(Function<Matrix, Vector> &fun, Vector &x, IterationState &state, int steps) {
            auto &f_gradient = state.f_gradient;
            auto &f_hessian = state.f_hessian;
            auto &f_diag_hessian = state.f_diag_hessian;

            auto &b_gradient = state.b_gradient;
            auto &b_diag_hessian = state.b_diag_hessian;

            auto &correction = state.correction;
            auto &diag_op = state.diag_op;
            auto &residual = state.residual;

            for (int ps = 0; ps < steps; ++ps) {
                if (use_non_linear_residual_) {
                    fun.gradient(x, f_gradient);
                }

                residual = f_gradient;

                if (!b_gradient.empty()) {
                    b_gradient.set(0.0);
                } else {
                    b_gradient.zeros(layout(x));
                }

                barrier_->gradient(x, b_gradient);

                residual += b_gradient;
                residual *= -1;

                barrier_->hessian_diag(x, b_diag_hessian);

                diag_op = f_diag_hessian + b_diag_hessian;
                diag_op = 1. / diag_op;

                correction = e_mul(diag_op, residual);

                handle_boundary(correction);

                Scalar alpha = line_search_projection_->compute(x, correction);

                state.alpha_stats.alpha_found(alpha);

                x += alpha * correction;

                if (!use_non_linear_residual_) {
                    f_gradient += alpha * (f_hessian * correction);
                }
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_BARRIER_MULTIGRID_HPP
