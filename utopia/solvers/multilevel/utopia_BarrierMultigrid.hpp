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

#include "utopia_Agglomerate.hpp"
#include "utopia_MatrixAgglomerator.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"

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
        using Traits = utopia::Traits<Vector>;
        using Layout = typename Traits::Layout;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;

        using Super = utopia::NonLinearSolver<Vector>;

        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;
        using Smoother = utopia::IterativeSolver<Matrix, Vector>;
        using SmootherPtr = std::shared_ptr<Smoother>;
        using IterativeSolver = utopia::IterativeSolver<Matrix, Vector>;
        using Level = utopia::Level<Matrix, Vector>;
        using Transfer = utopia::Transfer<Matrix, Vector>;
        using IPTransfer = utopia::IPTransfer<Matrix, Vector>;
        using VariableBoundSolverInterface = utopia::VariableBoundSolverInterface<Vector>;

        // using LogBarrierBase = utopia::LogBarrierBase<Matrix, Vector>;

        using LogBarrierFunctionBase = utopia::LogBarrierFunctionBase<Matrix, Vector>;

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

        class LevelMemory {
        public:
            std::shared_ptr<Matrix> matrix;
            Vector residual;
            Vector correction;
            Vector diag;
            Vector work;
            Vector modified_diag_op;

            Vector barrier_diag;
            bool matrix_changed{true};

            std::unique_ptr<Smoother> linear_smoother;

            void ensure() {
                if (!matrix) {
                    matrix = std::make_shared<Matrix>();
                }
            }
        };

        class NonlinearIterationState {
        public:
            // Function quantities
            Vector function_gradient;

            // Barrier quantities
            Vector barrier_gradient;

            // Constraints and boundary conditions
            Vector contraints_mask;

            // Stats
            AlphaStats alpha_stats;
        };

    public:
        BarrierMultigrid(const std::shared_ptr<LinearSolver> &coarse_solver) : coarse_solver_(coarse_solver) {}

        ~BarrierMultigrid() override = default;

        void read(Input &in) override {
            Super::read(in);

            // std::string barrier_function_type = "BoundedLogBarrier";
            std::string barrier_function_type = "BoundedLogBarrierFunction";

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
                .add_option("keep_initial_coarse_spaces",
                            keep_initial_coarse_spaces_,
                            "Compute the transfer operators only the first iteration.")
                .add_option("amg_n_coarse_spaces",
                            amg_n_coarse_spaces_,
                            "If an agglomerator is provide, define how many coarse spaces have to be generated,")
                .add_option("line_search_projection_weight",
                            non_smooth_projection_weight_,
                            "Coarse grid correction is projected in a non-smooth way. Value in [0, 1].")
                .parse(in);

            // barrier_ = LogBarrierFunctionFactory<Matrix, Vector>::new_log_barrier(barrier_function_type);
            barrier_ = LogBarrierFunctionFactory<Matrix, Vector>::new_log_barrier_function(barrier_function_type);
            barrier_->read(in);
        }

        void set_transfer_operators(const std::vector<std::shared_ptr<Transfer>> &transfer_operators) {
            this->transfer_operators_ = transfer_operators;
        }

        bool update_hessian(Function<Matrix, Vector> &fun,
                            const Vector &x,
                            NonlinearIterationState &state,
                            LevelMemory &mem) {
            if (!fun.hessian(x, *mem.matrix)) {
                return false;
            }

            if (empty(state.contraints_mask)) {
                generate_mask_from_matrix(*mem.matrix, state.contraints_mask, 1, 0);
            }

            if (agglomerator_) {
                init_algebraic();
            } else {
                init_with_user_transfer_operators();
            }

            return true;
        }

        bool solve(Function<Matrix, Vector> &fun, Vector &x) {
            bool converged = false;

            std::string mg_header_message = "BarrierMultigrid: " + std::to_string(n_levels()) + " levels";
            this->init_solver(mg_header_message, {" it. ", "|| g ||", "||g||/||g_0||"});

            ensure_defaults();

            NonlinearIterationState state;
            int top_level = n_levels() - 1;
            auto &mem = memory_[top_level];

            Scalar g_norm_0 = compute_norm_gradient_objective(fun, x, state, mem);
            PrintInfo::print_iter_status(0, {g_norm_0, 1});

            if (fun.is_hessian_constant()) {
                update_hessian(fun, x, state, mem);
            }

            for (SizeType it = 1; it < this->max_it() && !converged; ++it) {
                if (!fun.is_hessian_constant()) {
                    update_hessian(fun, x, state, mem);
                }

                //////////////////////////////////////////////////////////////////////

                if (!use_non_linear_residual_) {
                    fun.gradient(x, state.function_gradient);
                }

                //////////////////////////////////////////////////////////////////////

                for (int s = 0; s < mg_steps_; ++s) {
                    ////////////////////////////////////////////////////////////
                    // Pre-smoothing
                    ////////////////////////////////////////////////////////////
                    nonlinear_smooth(fun, x, state, pre_smoothing_steps());

                    ////////////////////////////////////////////////////////////
                    // Coarse grid correction
                    ////////////////////////////////////////////////////////////
                    if (use_coarse_space_) {
                        coarse_correction(fun, x, state);
                    }

                    ////////////////////////////////////////////////////////////
                    // Post-smoothing
                    ////////////////////////////////////////////////////////////
                    nonlinear_smooth(fun, x, state, post_smoothing_steps());
                }

                barrier_->update_barrier();

                /////////////////////////////////////////////
                // Convergence check

                Scalar g_norm = compute_norm_gradient_objective(fun, x, state, mem);
                PrintInfo::print_iter_status(it, {g_norm, g_norm / g_norm_0});
                converged = this->check_convergence(it, g_norm, (g_norm / g_norm_0), 1);

                /////////////////////////////////////////////
            }

            if (this->verbose()) {
                std::stringstream ss;
                state.alpha_stats.describe(ss);

                x.comm().root_print(ss.str());
            }

            if (debug_ && Traits::Backend == PETSC) {
                rename("r", state.barrier_gradient);

                write("R.m", state.barrier_gradient);

                rename("g", state.function_gradient);
                state.function_gradient *= -1;
                write("G.m", state.function_gradient);
            }

            clear_memory();
            return converged;
        }

        Scalar compute_norm_gradient_objective(Function<Matrix, Vector> &fun,
                                               const Vector &x,
                                               NonlinearIterationState &state,
                                               LevelMemory &mem) {
            fun.gradient(x, state.function_gradient);

            if (empty(state.barrier_gradient)) {
                state.barrier_gradient.zeros(layout(x));
            } else {
                state.barrier_gradient.set(0.);
            }

            barrier_->gradient(x, state.barrier_gradient);

            mem.residual = state.function_gradient + state.barrier_gradient;

            Scalar g_norm = norm2(mem.residual);
            return g_norm;
        }

        void init_memory(const Layout &) override { ensure_memory(); }

        int n_levels() const { return transfer_operators_.size() + 1; }

        void set_agglomerator(const std::shared_ptr<MatrixAgglomerator<Matrix>> &agglomerator) {
            agglomerator_ = agglomerator;
        }

        inline void set_barrier(const std::shared_ptr<LogBarrierFunctionBase> &barrier) { barrier_ = barrier; }
        inline const std::shared_ptr<LogBarrierFunctionBase> &barrier() { return barrier_; }

    public:
        BarrierMultigrid *clone() const override {
            auto ret = utopia::make_unique<BarrierMultigrid>(std::shared_ptr<LinearSolver>(coarse_solver_->clone()));
            if (linear_smoother_clonable_) {
                ret->linear_smoother_clonable_ = std::shared_ptr<Smoother>(linear_smoother_clonable_->clone());
            }

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

        void set_linear_smoother(const std::shared_ptr<Smoother> &smoother) { linear_smoother_clonable_ = smoother; }

        UTOPIA_NVCC_PRIVATE

        std::shared_ptr<LinearSolver> coarse_solver_;
        std::vector<std::shared_ptr<Transfer>> transfer_operators_;
        std::vector<LevelMemory> memory_;
        std::shared_ptr<LogBarrierFunctionBase> barrier_;
        std::shared_ptr<LineSearchBoxProjection<Vector>> line_search_projection_;
        std::shared_ptr<Smoother> linear_smoother_clonable_;

        bool debug_{true};
        bool use_coarse_space_{false};
        bool use_non_linear_residual_{false};
        bool keep_initial_coarse_spaces_{false};

        int mg_steps_{1};
        int pre_smoothing_steps_{3};
        int post_smoothing_steps_{3};

        int amg_n_coarse_spaces_{1};

        Scalar dumping_{0.98};
        Scalar non_smooth_projection_weight_{0.9};
        std::shared_ptr<MatrixAgglomerator<Matrix>> agglomerator_;

        void ensure_defaults() {
            if (!barrier_) {
                // barrier_ = std::make_shared<BoundedLogBarrier<Matrix, Vector>>();
                barrier_ = std::make_shared<BoundedLogBarrierFunction<Matrix, Vector>>();
            }

            if (!line_search_projection_) {
                line_search_projection_ = std::make_shared<LineSearchBoxProjection<Vector>>();
            }

            ensure_memory();

            auto &&box = this->get_box_constraints();
            barrier_->set_box_constraints(make_ref(box));
            barrier_->reset();
            line_search_projection_->set_box_constraints(make_ref(box));
        }

        void ensure_memory() {
            if (agglomerator_) {
                if (memory_.empty() || static_cast<int>(memory_.size()) != amg_n_coarse_spaces_ + 1) {
                    memory_.resize(amg_n_coarse_spaces_ + 1);
                    transfer_operators_.resize(amg_n_coarse_spaces_);
                }
            } else {
                if (memory_.size() != transfer_operators_.size() + 1) {
                    memory_.resize(transfer_operators_.size() + 1);
                }
            }

            for (auto &mem : memory_) {
                mem.ensure();

                if (linear_smoother_clonable_) {
                    mem.linear_smoother = std::unique_ptr<Smoother>(linear_smoother_clonable_->clone());
                }
            }
        }

        void clear_memory() { memory_.clear(); }

        void handle_eq_constraints(const NonlinearIterationState &state, Vector &v) const {
            v = e_mul(state.contraints_mask, v);
        }

        void linear_mg(const int l) {
            auto &&mem = memory_[l];
            auto &&transfer = transfer_operators_[l];

            transfer->restrict(memory_[l + 1].residual, mem.residual);
            transfer->restrict(memory_[l + 1].barrier_diag, mem.barrier_diag);

            if (empty(mem.correction)) {
                mem.correction.zeros(layout(mem.residual));
            } else {
                mem.correction.set(0.);
            }

            if (l == 0) {
                mem.matrix->shift_diag(mem.barrier_diag);
                coarse_solver_->solve(*mem.matrix, mem.residual, mem.correction);
                mem.matrix->set_diag(mem.diag);
            } else {
                update_linear_smoother_pre(mem);

                linear_smooth(pre_smoothing_steps(), mem);

                linear_mg(l - 1);

                linear_smooth(post_smoothing_steps(), mem);

                update_linear_smoother_post(mem);
            }

            transfer->interpolate(memory_[l].correction, memory_[l + 1].work);
            memory_[l + 1].correction += memory_[l + 1].work;
        }

        void update_linear_smoother_pre(LevelMemory &mem) {
            bool disable = true;
            if (!disable && mem.linear_smoother) {
                mem.matrix->shift_diag(mem.barrier_diag);
                mem.linear_smoother->update(mem.matrix);
            } else {
                // Jacobi
                mem.modified_diag_op = mem.diag + mem.barrier_diag;
                mem.modified_diag_op = 1. / mem.modified_diag_op;
            }
        }

        void update_linear_smoother_post(LevelMemory &mem) {
            if (mem.linear_smoother) {
                mem.matrix->set_diag(mem.diag);
            }
        }

        void linear_smooth(const int steps, LevelMemory &mem) {
            bool disable = true;
            if (!disable && mem.linear_smoother) {
                mem.linear_smoother->sweeps(steps);

                if (empty(mem.work)) {
                    mem.work.zeros(layout(mem.residual));
                } else {
                    mem.work.set(0.);
                }

                mem.linear_smoother->smooth(mem.residual, mem.work);

                // Update correction
                mem.correction += mem.work;

                // Update residual
                mem.residual -= *mem.matrix * mem.work;
                mem.residual -= e_mul(mem.barrier_diag, mem.work);

            } else {
                for (int ps = 0; ps < steps; ++ps) {
                    // Jacobi
                    mem.work = e_mul(mem.modified_diag_op, mem.residual);

                    // Update correction
                    mem.correction += mem.work;

                    // Update residual
                    mem.residual -= *mem.matrix * mem.work;
                    mem.residual -= e_mul(mem.barrier_diag, mem.work);
                }
            }
        }

        void coarse_correction(Function<Matrix, Vector> &fun, Vector &x, NonlinearIterationState &state) {
            ////////////////////////////////////////////////////////////////////////////////
            {
                // Refresh barrier function quantities

                int top_level = n_levels() - 1;
                auto &mem = memory_[top_level];

                mem.barrier_diag *= 0.;
                barrier_->hessian_diag(x, mem.barrier_diag);

                if (use_non_linear_residual_) {
                    fun.gradient(x, state.function_gradient);
                }

                mem.residual = state.function_gradient;
                barrier_->gradient(x, mem.residual);
                mem.residual *= -1;

                handle_eq_constraints(state, mem.residual);

                mem.correction.set(0.0);
            }

            ////////////////////////////////////////////////////////////////////////////////

            const int n_coarse_levels = this->n_levels() - 1;

            linear_mg(n_coarse_levels - 1);

            //////////////////////////////////////////////////////////////////////

            {
                int top_level = n_levels() - 1;
                auto &mem = memory_[top_level];

                if (non_smooth_projection_weight_ != 0) {
                    mem.work = x + mem.correction;
                    barrier_->project_onto_feasibile_region(mem.work);
                    mem.work -= x;
                    mem.work *= dumping_;
                } else {
                    mem.work.set(0.);
                }

                if (non_smooth_projection_weight_ != 1.0) {
                    Scalar alpha = line_search_projection_->compute(x, mem.correction);

                    state.alpha_stats.alpha_found(alpha);

                    mem.work *= non_smooth_projection_weight_;
                    mem.work += ((1 - non_smooth_projection_weight_) * alpha) * mem.correction;
                }

                x += mem.work;

                if (!use_non_linear_residual_) {
                    state.function_gradient += ((*mem.matrix) * mem.work);
                }
            }
        }

        void nonlinear_smooth(Function<Matrix, Vector> &fun, Vector &x, NonlinearIterationState &state, int steps) {
            int top_level = n_levels() - 1;

            auto &mem = memory_[top_level];

            for (int ps = 0; ps < steps; ++ps) {
                if (use_non_linear_residual_) {
                    fun.gradient(x, state.function_gradient);
                }

                mem.residual = state.function_gradient;

                if (!state.barrier_gradient.empty()) {
                    state.barrier_gradient.set(0.0);
                } else {
                    state.barrier_gradient.zeros(layout(x));
                }

                barrier_->gradient(x, state.barrier_gradient);

                mem.residual += state.barrier_gradient;
                mem.residual *= -1;

                barrier_->hessian_diag(x, mem.barrier_diag);

                handle_eq_constraints(state, mem.residual);
                handle_eq_constraints(state, mem.barrier_diag);

                //////////////////////////////////////////////////
                if (mem.linear_smoother) {
                    mem.linear_smoother->sweeps(1);
                    mem.matrix->shift_diag(mem.barrier_diag);
                    mem.linear_smoother->update(mem.matrix);

                    if (empty(mem.correction)) {
                        mem.correction.zeros(layout(mem.residual));
                    } else {
                        mem.correction.set(0.);
                    }

                    mem.linear_smoother->smooth(mem.residual, mem.correction);

                    mem.matrix->set_diag(mem.diag);
                } else {
                    // Jacobi
                    mem.work = mem.diag + mem.barrier_diag;
                    mem.work = 1. / mem.work;

                    mem.correction = e_mul(mem.work, mem.residual);
                }

                //////////////////////////////////////////////////

                Scalar alpha = line_search_projection_->compute(x, mem.correction);

                state.alpha_stats.alpha_found(alpha);

                x += alpha * mem.correction;

                if (!use_non_linear_residual_) {
                    state.function_gradient += alpha * (*mem.matrix * mem.correction);
                }
            }
        }

        static void generate_mask_from_matrix(const Matrix &A,
                                              Vector &mask,
                                              const Scalar on_value,
                                              const Scalar off_value) {
            const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon() * 1e6;

            // FIXME once atomic_store is avaialable
            mask.values(row_layout(A), off_value);

            {
                auto mask_view = view_device(mask);

                A.read(UTOPIA_LAMBDA(const SizeType &i, const SizeType &j, const Scalar &value) {
                    if (i == j) return;

                    if (device::abs(value) > off_diag_tol) {
                        // mask_view.atomic_add(i, 1.0);
                        mask_view.set(i, on_value);
                    }
                });
            }
        }

        void init_with_user_transfer_operators() {
            int n_coarse_levels = n_levels() - 1;

            for (int i = n_coarse_levels - 1; i >= 0; --i) {
                auto &mem = memory_[i];
                transfer_operators_[i]->restrict(*memory_[i + 1].matrix, *mem.matrix);
            }

            for (auto &mem : memory_) {
                mem.diag = diag(*mem.matrix);
                mem.matrix_changed = true;
            }
        }

        void init_algebraic() {
            UTOPIA_TRACE_REGION_BEGIN("BarrierMultigrid::init_algebraic");

            if (amg_n_coarse_spaces_ == 0) {
                Utopia::Abort("BarrierMultigrid: Invalid number of coarse levels!");
            }

            //////////////////////////////////////////////////////////////////

            if (!keep_initial_coarse_spaces_ || !transfer_operators_[0]) {
                auto fine_matrix = memory_[amg_n_coarse_spaces_].matrix;

                for (SizeType l = amg_n_coarse_spaces_ - 1; l >= 0; --l) {
                    auto &mem = memory_[l];

                    assert(bool(fine_matrix));
                    assert(!empty(*fine_matrix));
                    assert(bool(mem.matrix));

                    auto t = agglomerator_->create_transfer(*fine_matrix);
                    t->restrict(*fine_matrix, *mem.matrix);

                    // disp(static_cast<MatrixTransfer<Matrix, Vector> &>(*t).I());

                    transfer_operators_[l] = t;
                    fine_matrix = mem.matrix;
                }
            } else if (keep_initial_coarse_spaces_) {
                auto fine_matrix = memory_[amg_n_coarse_spaces_].matrix;

                for (SizeType l = amg_n_coarse_spaces_ - 1; l >= 0; --l) {
                    auto &mem = memory_[l];

                    transfer_operators_[l]->restrict(*fine_matrix, *mem.matrix);

                    fine_matrix = mem.matrix;
                }
            }

            //////////////////////////////////////////////////////////////////

            for (auto &mem : memory_) {
                mem.diag = diag(*mem.matrix);
                mem.matrix_changed = true;
            }

            UTOPIA_TRACE_REGION_END("BarrierMultigrid::init_algebraic");
        }
    };  // namespace utopia

}  // namespace utopia

#endif  // UTOPIA_BARRIER_MULTIGRID_HPP
