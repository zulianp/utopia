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
        }

        void set_transfer_operators(const std::vector<std::shared_ptr<Transfer>> &transfer_operators) {
            auto last = transfer_operators.back();

            auto tt = std::dynamic_pointer_cast<TruncatedTransfer>(last);

            if (!tt) {
                Utopia::Abort(
                    "MonotoneMultigrid::set_transfer_operators(...) last operators must be of type TruncatedTransfer");
            }

            truncated_transfer_ = tt;
            transfer_operators.resize(transfer_operators.size() - 1);

            algo_.set_transfer_operators(transfer_operators);
        }

        void set_transfer_operators(const std::vector<std::shared_ptr<Transfer>> &transfer_operators,
                                    const std::shared_ptr<TruncatedTransfer> &last_level) {
            algo_.set_transfer_operators(transfer_operators);
            truncated_transfer_ = last_level;
        }

        void init(const std::shared_ptr<QPSmoother> &qp_smoother) { qp_smoother_ = fine_smoother; }

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

            if (algo_.must_perform_galerkin_assembly()) {
                std::shared_ptr<Matrix> op_H = std::make_shared<Matrix>();
                truncated_transfer_->restrict(*op, *op_H);

                if (algo_.fix_semidefinite_operators()) {
                    algo_.fix_semidefinite_operator(*op_H);
                }

                algo_.update(op_H);
            }
        }

        void init_memory(const Layout &l) override {
            Super::init_memory(l);
            active_set_->init(l);
        }

        inline bool step(const Vector &rhs, Vector &x, const SizeType &smoothing_steps = 1) {
            assert(l != this->n_levels() - 2);

            qp_smoother_->sweeps(smoothing_steps);

            if (auto *fine_smoother_vi = dynamic_cast<VariableBoundSolverInterface *>(fine_smoother)) {
                fine_smoother_vi->set_box_constraints(this->get_box_constraints());
            } else {
                utopia_error(
                    "MonotoneMultigrid: requires VariableBoundSolverInterface to be the fine level "
                    "smoother ");
            }

            fine_smoother->smooth(rhs, x);

            if (pre_sm) {
                if (auto *trunc_transfer =
                        dynamic_cast<IPTruncatedTransfer<Matrix, Vector> *>(this->transfers_[l - 1].get())) {
                    ////////////////////////////////////////////////////////////////
                    // FIXME duplicated code (ProjectedGaussSeidelQR needs changes to avoid this)
                    if (auto *pgs_QR = dynamic_cast<ProjectedGaussSeidelQR<Matrix, Vector> *>(fine_smoother)) {
                        trunc_transfer->truncate_interpolation(pgs_QR->get_active_set());

                        if (this->verbose()) {
                            Scalar n_active = sum(pgs_QR->get_active_set());
                            x.comm().root_print("n_active: " + std::to_string(SizeType(n_active)));
                        }

                        this->galerkin_assembly(this->get_operator());

                        for (std::size_t l = 1; l != smoothers_.size() - 1; ++l) {
                            smoothers_[l]->update(level(l).A_ptr());
                        }

                        coarse_solver_->update(level(0).A_ptr());
                        ////////////////////////////////////////////////////////////////
                    } else {
                        active_set_->verbose(this->verbose());
                        if (active_set_->determine(this->get_box_constraints(), x)) {
                            trunc_transfer->truncate_interpolation(active_set_->indicator());

                            this->galerkin_assembly(this->get_operator());

                            for (std::size_t l = 1; l != smoothers_.size() - 1; ++l) {
                                smoothers_[l]->update(level(l).A_ptr());
                            }

                            coarse_solver_->update(level(0).A_ptr());
                        }
                    }
                }

                else {
                    Utopia::Abort("MonotoneMultigrid: requires IPTruncatedTransfer for the finest level.");
                }
            }
            return true;
        }

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
    };

}  // namespace utopia

#endif  // UTOPIA_MONOTONE_MULTIGRID_HPP
