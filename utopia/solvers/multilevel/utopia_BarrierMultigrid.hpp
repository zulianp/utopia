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

#include "utopia_ActiveSet.hpp"
#include "utopia_Multigrid.hpp"

#include "utopia_LogBarrierFunctionBase.hpp"

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
    class BarrierMultigrid final : public QPSolver<Matrix, Vector> {
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

        class NonLinearSmoother : public Configurable {
            void update(const std::shared_ptr<const Matrix> &op) {}
            bool apply(const Vector &rhs, Vector &x) { return false; }

            void read(Input &in) override {}
        };

        using DefaultSmoother = NonLinearSmoother;

    public:
        BarrierMultigrid(const std::shared_ptr<LinearSolver> &coarse_solver,
                         const std::shared_ptr<NonLinearSmoother> &nl_smoother = std::make_shared<DefaultSmoother>())
            : coarse_solver_(coarse_solver), nl_smoother_(nl_smoother) {}

        ~BarrierMultigrid() override = default;

        void read(Input &in) override {
            Super::read(in);

            if (nl_smoother_) {
                in.get("nl_smoother", *nl_smoother_);
            }

            in.get("debug", debug_);
            in.get("use_coarse_space", use_coarse_space_);
            in.get("pre_smoothing_steps", pre_smoothing_steps_);
            in.get("post_smoothing_steps", post_smoothing_steps_);
        }

        void set_transfer_operators(const std::vector<std::shared_ptr<Transfer>> &transfer_operators) {
            this->transfer_operators_ = transfer_operators;
        }

        void init(const std::shared_ptr<NonLinearSmoother> &nl_smoother) { nl_smoother_ = nl_smoother; }

        void update(const std::shared_ptr<const Matrix> &op) override { Super::update(op); }

        bool apply(const Vector &rhs, Vector &x) override {
            bool converged = false;

            std::string mg_header_message = "BarrierMultigrid: " + std::to_string(n_levels()) + " levels";
            this->init_solver(mg_header_message, {" it. ", "|| P(x-g)-x ||", "|| u - u_old ||", "# active sets"});

            Vector x_old = x;
            Vector grad;
            ///
            grad = grad - rhs;
            Scalar g_norm = this->criticality_measure_inf(x, grad);

            if (this->verbose()) {
                PrintInfo::print_iter_status(0, {g_norm, 0});
            }

            for (SizeType it = 1; it < this->max_it() && !converged; ++it) {
                if (!step(rhs, x)) {
                    Utopia::Abort("BarrierMultigrid::step FAILED");
                }

                converged = this->check_convergence(it, g_norm, 1, 1);

                if (!converged) {
                    x_old = x;
                }
            }

            return converged;
        }

        void init_memory(const Layout &l) override { Super::init_memory(l); }

        inline bool step(const Vector &rhs, Vector &x) {
            // nl_smoother_->set_box_constraints(this->get_box_constraints());
            // nl_smoother_->sweeps(pre_smoothing_steps());
            // nl_smoother_->smooth(rhs, x);

            // TODO

            // nl_smoother_->sweeps(post_smoothing_steps());
            // nl_smoother_->smooth(rhs, x);
            return true;
        }

        void pre_smoothing_steps(const SizeType n) { pre_smoothing_steps(n); }
        void post_smoothing_steps(const SizeType n) { post_smoothing_steps(n); }

        int n_levels() const { return 2; }

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
            help_ = x - g;

            this->project(*this->constraints_.lower_bound(), *this->constraints_.upper_bound(), help_);
            help_ -= x;

            return norm2(help_);
        }

        inline int pre_smoothing_steps() const { return pre_smoothing_steps_; }
        inline int post_smoothing_steps() const { return post_smoothing_steps_; }

    protected:
        std::shared_ptr<LinearSolver> coarse_solver_;
        std::shared_ptr<NonLinearSmoother> nl_smoother_;

        bool debug_{false};
        bool use_coarse_space_{true};

        int pre_smoothing_steps_{3};
        int post_smoothing_steps_{3};

        Vector help_;

    };  // namespace utopia

}  // namespace utopia

#endif  // UTOPIA_BARRIER_MULTIGRID_HPP
