#ifndef UTOPIA_MONOTONE_ALGEBRAIC_MULTIGRID_HPP
#define UTOPIA_MONOTONE_ALGEBRAIC_MULTIGRID_HPP

#include "utopia_Agglomerate.hpp"
#include "utopia_AlgebraicMultigrid.hpp"
#include "utopia_Clonable.hpp"
#include "utopia_IPRTransfer.hpp"
#include "utopia_MonotoneMultigrid.hpp"
#include "utopia_QPSolver.hpp"

namespace utopia {

    template <class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class MonotoneAlgebraicMultigrid final : public QPSolver<Matrix, Vector> {
    public:
        using Solver = utopia::LinearSolver<Matrix, Vector>;
        using Smoother = utopia::IterativeSolver<Matrix, Vector>;
        using SmootherPtr = std::shared_ptr<Smoother>;
        using IterativeSolver = utopia::IterativeSolver<Matrix, Vector>;
        using Level = utopia::Level<Matrix, Vector>;
        using Transfer = utopia::Transfer<Matrix, Vector>;
        using VariableBoundSolverInterface = utopia::VariableBoundSolverInterface<Vector>;
        using Super = utopia::QPSolver<Matrix, Vector>;
        using Layout = typename Traits<Vector>::Layout;
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Agglomerator = utopia::MatrixAgglomerator<Matrix>;
        using QPSmoother = utopia::QPSolver<Matrix, Vector>;

        using Super::max_it;
        using Super::verbose;

        void init() { algorithm_.algorithm().fix_semidefinite_operators(true); }

        MonotoneAlgebraicMultigrid()
            : qp_smoother_(std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>()),
              active_set_(std::make_shared<ActiveSet<Vector>>()),
              agglomerator_(std::make_shared<Agglomerate<Matrix>>()),
              algorithm_(std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>(),
                         std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>(),
                         agglomerator_) {
            init();
        }

        MonotoneAlgebraicMultigrid(const std::shared_ptr<QPSmoother> &qp_smoother,
                                   const std::shared_ptr<Smoother> &smoother,
                                   const std::shared_ptr<Solver> &coarse_solver,
                                   const std::shared_ptr<MatrixAgglomerator<Matrix>> &agglomerator)
            : qp_smoother_(qp_smoother),
              active_set_(std::make_shared<ActiveSet<Vector>>()),
              agglomerator_(agglomerator),
              algorithm_(smoother, coarse_solver, agglomerator) {
            init();
        }

        MonotoneAlgebraicMultigrid(const std::shared_ptr<QPSmoother> &qp_smoother,
                                   const std::shared_ptr<Solver> &coarse_solver)
            : qp_smoother_(qp_smoother),
              active_set_(std::make_shared<ActiveSet<Vector>>()),
              agglomerator_(std::make_shared<Agglomerate<Matrix>>()),
              algorithm_(qp_smoother->clone(), coarse_solver, agglomerator_) {
            init();
        }

        ~MonotoneAlgebraicMultigrid() override = default;

        inline static std::shared_ptr<IPRTruncatedTransfer<Matrix, Vector>> new_fine_level_transfer(
            const std::shared_ptr<Matrix> &interpolation) {
            return std::make_shared<IPRTruncatedTransfer<Matrix, Vector>>(interpolation);
        }

        inline static std::shared_ptr<IPRTransfer<Matrix, Vector>> new_coarse_level_transfer(
            const std::shared_ptr<Matrix> &interpolation) {
            return std::make_shared<IPRTransfer<Matrix, Vector>>(interpolation);
        }

        void verbose(const bool &val) override {
            Super::verbose(val);
            algorithm_.verbose(val);
        }

        void max_it(const SizeType &max_it_in) override {
            Super::max_it(max_it_in);
            algorithm_.max_it(max_it_in);
        }

        void read(Input &in) override {
            Super::read(in);
            algorithm_.read(in);
            in.get("active_set", *active_set_);
        }

        void update(const std::shared_ptr<const Matrix> &op) override {
            Super::update(op);
            init_memory(row_layout(*op));
            modified_op_ = nullptr;
            qp_smoother_->update(op);
        }

        void init_memory(const Layout &l) override {
            Super::init_memory(l);
            active_set_->init(l);
        }

        inline void set_n_levels(const int n_levels) { algorithm_.set_n_levels(std::max(2, n_levels - 1)); }

        /**
         * @brief      The solve function for multigrid method.
         *
         * @param[in]  rhs   The right hand side.
         * @param      x   The initial guess.
         *
         */
        bool apply(const Vector &rhs, Vector &x) override {
            if (this->verbose()) {
                this->init_solver("utopia::MonotoneAlgebraicMultigrid", {" it. ", "|| u - u_old ||"});
            }

            Vector x_old = x;
            Vector diff;

            bool converged = false;
            for (int it = 0; it < this->max_it() && !converged; ++it) {
                apply_once(rhs, x);

                diff = x_old - x;
                Scalar diff_norm = norm2(diff);
                x_old = x;

                if (this->verbose()) {
                    PrintInfo::print_iter_status(it, {diff_norm});
                }

                converged = this->check_convergence(it, 1, 1, diff_norm);
            }

            return converged;
        }

        bool apply_once(const Vector &rhs, Vector &x) {
            auto &op = *this->get_operator();

            bool first = false;
            if (!modified_op_) {
                // Make a copy
                modified_op_ = std::make_shared<Matrix>(op);
                first = true;
            }

            auto &&box = this->get_box_constraints();

            qp_smoother_->set_box_constraints(box);

            // Pre-smoothing
            qp_smoother_->smooth(rhs, x);

            if (active_set_->determine(box, x)) {
                if (!first) {
                    modified_op_->same_nnz_pattern_copy(op);
                }

                set_zero_rows(*modified_op_, active_set_->indicator(), 1.0);
                transfer_ = agglomerator_->create_transfer(*modified_op_);
                auto coarse_mat = std::make_shared<Matrix>();
                transfer_->restrict(*modified_op_, *coarse_mat);
                algorithm_.algorithm().fix_semidefinite_operator(*coarse_mat);
                // rename("a_H", *coarse_mat);
                // write("A_H.m", *coarse_mat);

                algorithm_.update(coarse_mat);

            } else {
                if (first) {
                    transfer_ = agglomerator_->create_transfer(*modified_op_);
                    auto coarse_mat = std::make_shared<Matrix>();
                    transfer_->restrict(*modified_op_, *coarse_mat);
                    algorithm_.update(coarse_mat);
                }
            }

            residual = rhs - (*modified_op_) * x;
            active_set_->zero_out_active(residual);

            transfer_->restrict(residual, coarse_residual);

            coarse_correction.zeros(layout(coarse_residual));

            algorithm_.smooth(coarse_residual, coarse_correction);

            transfer_->interpolate(coarse_correction, correction);

            // See if necessary
            active_set_->zero_out_active(correction);

            x += correction;

            // Post-smoothing
            qp_smoother_->smooth(rhs, x);

            return true;
        }

        MonotoneAlgebraicMultigrid *clone() const override {
            auto amg = utopia::make_unique<MonotoneAlgebraicMultigrid>();
            amg->qp_smoother_ = std::shared_ptr<QPSmoother>(qp_smoother_->clone());
            // amg->active_set_ = std::shared_ptr<ActiveSet<Vector>>(active_set_->clone());
            amg->agglomerator_ = std::shared_ptr<Agglomerator>(agglomerator_->clone());
            return amg.release();
        }

        inline ActiveSet<Vector> &active_set() { return *active_set_; }

        inline void active_set(const std::shared_ptr<ActiveSet<Vector>> active_set) {
            active_set_ = std::move(active_set);
        }

    private:
        std::shared_ptr<QPSmoother> qp_smoother_;
        std::shared_ptr<ActiveSet<Vector>> active_set_;
        std::shared_ptr<Agglomerator> agglomerator_;
        AlgebraicMultigrid<Matrix, Vector> algorithm_;

        std::shared_ptr<Matrix> modified_op_;
        std::shared_ptr<Transfer> transfer_;

        Vector residual, correction, coarse_residual, coarse_correction;
    };

}  // namespace utopia

#endif  // UTOPIA_MONOTONE_ALGEBRAIC_MULTIGRID_HPP
