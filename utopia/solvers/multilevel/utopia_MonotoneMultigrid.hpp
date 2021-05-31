#ifndef UTOPIA_MONOTONE_MULTIGRID_HPP
#define UTOPIA_MONOTONE_MULTIGRID_HPP
#include "utopia_ConvergenceReason.hpp"
#include "utopia_Core.hpp"

#include "utopia_IPRTruncatedTransfer.hpp"
#include "utopia_IPTruncatedTransfer.hpp"
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
#include "utopia_ProjectedGaussSeidelNew.hpp"
#include "utopia_ProjectedGaussSeidelQR.hpp"

#include <cassert>
#include <ctime>

namespace utopia {
    /**
     * @brief      The class for  MonotoneMultigrid solver.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class MonotoneMultigrid final : public QPSolver<Matrix, Vector> {
        using Layout = typename Traits<Vector>::Layout;
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        using QPSmoother = utopia::QPSolver<Matrix, Vector>;
        using Super = utopia::QPSolver<Matrix, Vector>;

        using Solver = utopia::LinearSolver<Matrix, Vector>;
        using Smoother = utopia::IterativeSolver<Matrix, Vector>;
        using SmootherPtr = std::shared_ptr<Smoother>;
        using IterativeSolver = utopia::IterativeSolver<Matrix, Vector>;
        using Level = utopia::Level<Matrix, Vector>;
        using Transfer = utopia::Transfer<Matrix, Vector>;
        using TruncatedTransfer = utopia::IPTruncatedTransfer<Matrix, Vector>;

        using VariableBoundSolverInterface = utopia::VariableBoundSolverInterface<Vector>;
        using DefaultSmoother = utopia::ProjectedGaussSeidel<Matrix, Vector>;

    public:
        static const int V_CYCLE = 1;
        static const int W_CYCLE = 2;

        MonotoneMultigrid(const std::shared_ptr<QPSmoother> &qp_smoother,
                          const std::shared_ptr<Smoother> &coarse_smoother,
                          const std::shared_ptr<Solver> &coarse_solver)
            : algo_(coarse_smoother, coarse_solver),
              qp_smoother_(qp_smoother),
              active_set_(std::make_shared<ActiveSet<Vector>>()) {
            algo_.must_generate_masks(false);
            algo_.fix_semidefinite_operators(true);
        }

        MonotoneMultigrid()
            : algo_(std::make_shared<DefaultSmoother>(), std::make_shared<DefaultSmoother>()),
              qp_smoother_(std::make_shared<DefaultSmoother>()),
              active_set_(std::make_shared<ActiveSet<Vector>>()) {
            algo_.must_generate_masks(false);
            algo_.fix_semidefinite_operators(true);
        }

        ~MonotoneMultigrid() override = default;

        void read(Input &in) override {
            Super::read(in);
            algo_.read(in);

            if (active_set_) {
                in.get("active_set", *active_set_);
            }

            if (qp_smoother_) {
                in.get("qp_smoother", *qp_smoother_);
            }

            in.get("coarse_space_it", coarse_space_it_);
            in.get("debug", debug_);
            in.get("use_coarse_space", use_coarse_space_);

            algo_.verbose(false);
        }

        void set_transfer_operators(const std::vector<std::shared_ptr<Transfer>> &transfer_operators) {
            auto last = transfer_operators.back();

            auto tt = std::dynamic_pointer_cast<TruncatedTransfer>(last);

            if (!tt) {
                Utopia::Abort(
                    "MonotoneMultigrid::set_transfer_operators(...) last operators must be of type TruncatedTransfer");
            }

            truncated_transfer_ = tt;

            auto transfer_operators_copy = transfer_operators;
            transfer_operators_copy.resize(transfer_operators.size() - 1);

            algo_.set_transfer_operators(transfer_operators_copy);
        }

        void set_transfer_operators(const std::vector<std::shared_ptr<Transfer>> &transfer_operators,
                                    const std::shared_ptr<TruncatedTransfer> &last_level) {
            algo_.set_transfer_operators(transfer_operators);
            truncated_transfer_ = last_level;
        }

        void init(const std::shared_ptr<QPSmoother> &qp_smoother) { qp_smoother_ = qp_smoother; }

        inline static std::shared_ptr<TruncatedTransfer> new_fine_level_transfer(
            const std::shared_ptr<Matrix> &interpolation) {
            return std::make_shared<IPTruncatedTransfer<Matrix, Vector>>(interpolation);
        }

        inline static std::shared_ptr<TruncatedTransfer> new_coarse_level_transfer(
            const std::shared_ptr<Matrix> &interpolation) {
            return std::make_shared<IPTransfer<Matrix, Vector>>(interpolation);
        }

        /*! @brief if overriden the subclass has to also call this one first
         */
        void update(const std::shared_ptr<const Matrix> &op) override {
            Super::update(op);
            init_memory(row_layout(*op));
            qp_smoother_->update(op);

            if (algo_.must_perform_galerkin_assembly()) {
                galerkin_assembly(op);
            }
        }

        inline SizeType n_levels() const { return algo_.n_levels() + 1; }

        bool apply(const Vector &rhs, Vector &x) override {
            bool converged = false;

            // algo_.verbose(this->verbose());

            std::string mg_header_message = "MontoneMultigrid: " + std::to_string(n_levels()) + " levels";
            this->init_solver(mg_header_message, {" it. ", "|| u - u_old ||"});

            Vector x_old = x;
            for (SizeType it = 0; it < this->max_it() && !converged; ++it) {
                if (!step(rhs, x)) {
                    Utopia::Abort("MonotoneMultigrid::step FAILED");
                }

                Scalar x_diff = norm2(x - x_old);

                if (this->verbose()) {
                    PrintInfo::print_iter_status(it, {x_diff});
                }

                converged = this->check_convergence(it, 1, 1, x_diff);

                if (!converged) {
                    x_old = x;
                }
            }

            return converged;
        }

        void init_memory(const Layout &l) override {
            Super::init_memory(l);
            active_set_->init(l);
        }

        void galerkin_assembly(const std::shared_ptr<const Matrix> &op) {
            if (use_coarse_space_) {
                std::shared_ptr<Matrix> op_H = std::make_shared<Matrix>();
                truncated_transfer_->restrict(*op, *op_H);

                if (algo_.fix_semidefinite_operators()) {
                    algo_.fix_semidefinite_operator(*op_H);
                }

                if (n_levels() > 2) {
                    algo_.update(op_H);
                } else {
                    algo_.coarse_solver()->update(op_H);
                }
            }
        }

        void coarse_space_step(const Vector &r_coarse, Vector &c_coarse) {
            if (n_levels() > 2) {
                algo_.max_it(coarse_space_it_);
                algo_.apply(r_coarse, c_coarse);
            } else {
                algo_.coarse_solver()->apply(r_coarse, c_coarse);
            }
        }

        inline bool step(const Vector &rhs, Vector &x) {
            qp_smoother_->set_box_constraints(this->get_box_constraints());
            qp_smoother_->sweeps(algo_.pre_smoothing_steps());
            qp_smoother_->smooth(rhs, x);

            if (use_coarse_space_) {
                if (auto *pgs_QR = dynamic_cast<ProjectedGaussSeidelQR<Matrix, Vector> *>(qp_smoother_.get())) {
                    truncated_transfer_->truncate_interpolation(pgs_QR->get_active_set());

                    if (this->verbose()) {
                        Scalar n_active = sum(pgs_QR->get_active_set());
                        x.comm().root_print("n_active: " + std::to_string(SizeType(n_active)));
                    }

                    this->galerkin_assembly(this->get_operator());
                    ////////////////////////////////////////////////////////////////
                } else {
                    active_set_->verbose(this->verbose());

                    if (active_set_->determine(this->get_box_constraints(), x)) {
                        truncated_transfer_->truncate_interpolation(active_set_->indicator());
                        this->galerkin_assembly(this->get_operator());
                    }
                }

                Vector r, c;
                Vector r_coarse, c_coarse;

                r = (*this->get_operator()) * x;
                r = rhs - r;

                active_set_->zero_out_active(r);

                truncated_transfer_->restrict(r, r_coarse);

                if (debug_) {
                    Scalar r_norm = norm2(r);
                    Scalar r_coarse_norm = norm2(r_coarse);

                    std::stringstream ss;

                    ss << "residual before correction: " << r_norm << '\n';
                    ss << "resricted residual: " << r_coarse_norm << '\n';

                    x.comm().root_print(ss.str());
                }

                c_coarse.zeros(layout(r_coarse));

                // Compute coarse grid correction
                coarse_space_step(r_coarse, c_coarse);

                truncated_transfer_->interpolate(c_coarse, c);

                active_set_->zero_out_active(c);

                x += c;

                if (debug_) {
                    r = (*this->get_operator()) * x;
                    r = rhs - r;

                    Scalar r_norm = norm2(r);
                    Scalar c_norm = norm2(c);
                    Scalar c_coarse_norm = norm2(c_coarse);

                    std::stringstream ss;

                    ss << "residual after correction: " << r_norm << '\n';
                    ss << "norm correction: " << c_norm << '\n';
                    ss << "norm coarse correction: " << c_coarse_norm << '\n';

                    x.comm().root_print(ss.str());
                }
            }

            qp_smoother_->sweeps(algo_.post_smoothing_steps());
            qp_smoother_->smooth(rhs, x);
            return true;
        }

        void pre_smoothing_steps(const SizeType n) { algo_.pre_smoothing_steps(n); }
        void post_smoothing_steps(const SizeType n) { algo_.post_smoothing_steps(n); }

    public:
        MonotoneMultigrid *clone() const override {
            auto ret = utopia::make_unique<MonotoneMultigrid>();
            ret->algo_ = algo_;
            return ret.release();
        }

        /**
         * @brief      Function changes direct solver needed for coarse grid solve.
         *
         * @param[in]  linear_solver  The linear solver.
         *
         * @return
         */
        bool change_coarse_solver(const std::shared_ptr<Solver> &coarse_solver) {
            return algo_.change_coarse_solver(coarse_solver);
        }

        inline ActiveSet<Vector> &active_set() { return *active_set_; }

        inline void active_set(const std::shared_ptr<ActiveSet<Vector>> active_set) {
            active_set_ = std::move(active_set);
        }

    protected:
        Multigrid<Matrix, Vector> algo_;
        std::shared_ptr<QPSmoother> qp_smoother_;
        std::shared_ptr<TruncatedTransfer> truncated_transfer_;
        std::shared_ptr<ActiveSet<Vector>> active_set_;
        SizeType coarse_space_it_{1};
        bool debug_{false};
        bool use_coarse_space_{true};
    };  // namespace utopia

}  // namespace utopia

#endif  // UTOPIA_MONOTONE_MULTIGRID_HPP
