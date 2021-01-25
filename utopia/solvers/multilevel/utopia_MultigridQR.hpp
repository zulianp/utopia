#ifndef UTOPIA_MULTIGRID_QR_HPP
#define UTOPIA_MULTIGRID_QR_HPP
#include "utopia_ConvergenceReason.hpp"
#include "utopia_Core.hpp"
#include "utopia_IPRTruncatedTransfer.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Level.hpp"
#include "utopia_LinearMultiLevel.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_ProjectedGaussSeidelQR.hpp"
#include "utopia_Recorder.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_Utils.hpp"

#include <cassert>
#include <ctime>

namespace utopia {
    /**
     * @brief      The class for Linear Multigrid solver.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class MultigridQR final : public LinearMultiLevel<Matrix, Vector>, public VariableBoundSolverInterface<Vector> {
        typedef utopia::LinearSolver<Matrix, Vector> Solver;
        typedef utopia::IterativeSolver<Matrix, Vector> Smoother;
        using SmootherPtr = std::shared_ptr<Smoother>;

        typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
        typedef utopia::Level<Matrix, Vector> Level;
        typedef utopia::Transfer<Matrix, Vector> Transfer;
        using VariableBoundSolverInterface = utopia::VariableBoundSolverInterface<Vector>;

        using Super = utopia::LinearMultiLevel<Matrix, Vector>;

        using Layout = typename Traits<Vector>::Layout;
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        typedef struct {
            std::vector<Vector> r, x, rhs;

            void init(const int n_levels) {
                r.resize(n_levels);
                x.resize(n_levels);
                rhs.resize(n_levels);
            }

            inline bool valid(const std::size_t n_levels) const { return r.size() == n_levels; }

        } LevelMemoryMGQR;

        LevelMemoryMGQR memory;

        // #define BENCHMARKING_mode
        // #define CHECK_NUM_PRECISION_mode

    public:
        static const int V_CYCLE = 1;
        static const int W_CYCLE = 2;

        /**
         * @brief      Multigrid class
         *
         * @param[in]  smoother       The smoother.
         * @param[in]  coarse_solver  The direct solver for coarse level.
         */
        MultigridQR(const std::shared_ptr<Smoother> &fine_smoother,
                    const std::shared_ptr<Smoother> &coarse_smoother,
                    const std::shared_ptr<Solver> &coarse_solver,
                    const SizeType &num_levels)
            : coarse_solver_(coarse_solver), use_line_search_(false) {
            this->must_generate_masks(true);
            this->num_levels_ = num_levels;

            smoothers_.resize(this->n_levels());
            smoothers_[this->num_levels_ - 1] = fine_smoother;

            for (SizeType l = 1; l < this->num_levels_ - 1; ++l) {
                smoothers_[l] = std::shared_ptr<Smoother>(coarse_smoother->clone());
            }
        }

        ~MultigridQR() override = default;

        void read(Input &in) override {
            Super::read(in);

            in.get("use_line_search", use_line_search_);

            // TODO:: read params for all smoothers
            // if(smoother_cloneable_) {
            //     in.get("smoother", *smoother_cloneable_);
            // }

            if (coarse_solver_) {
                in.get("coarse_solver", *coarse_solver_);
            }
        }

        void print_usage(std::ostream &os) const override {
            Super::print_usage(os);

            this->print_param_usage(
                os, "perform_galerkin_assembly", "bool", "Flag turning on/off galerkin assembly.", "true");
            this->print_param_usage(os,
                                    "use_line_search",
                                    "bool",
                                    "Flag turning on/off line-search after coarse grid correction.",
                                    "false");
            this->print_param_usage(os, "block_size", "int", "Block size for systems.", "1");
            // this->print_param_usage(os, "smoother", "Smoother", "Input parameters for
            // all smoothers.", "-");
            this->print_param_usage(os, "coarse_solver", "LinearSolver", "Input parameters for coarse solver.", "-");
        }

        /*! @brief if overriden the subclass has to also call this one first
         */
        void update(const std::shared_ptr<const Matrix> &op) override {
            Super::update(op);
            this->galerkin_assembly(op);

            update();
        }

        /*! @brief if no galerkin assembly is used but instead set_linear_operators is
           used. One can call this update instead of the other one.
         */
        void update() {
            smoothers_.resize(this->n_levels());
            // smoothers_[0] = nullptr;

            for (std::size_t l = 1; l != smoothers_.size(); ++l) {
                if (smoothers_[l] == nullptr) {
                    assert(smoothers_[l]);
                }
                smoothers_[l]->update(level(l).A_ptr());
            }
            coarse_solver_->update(level(0).A_ptr());
        }

        // rememebr to put level -1, as c++ indexing goes from 0
        void set_smoother(const std::shared_ptr<Smoother> &smoother, const SizeType &level) {
            if (level > this->n_levels() || level < 1) {
                utopia_error("MultigridQR:: set_smoother, invalid level.");
            }

            if (smoothers_.size() - 1 < level) {
                utopia_error("MultigridQR:: set_smoother, smoothers array was not allocated.");
            }

            smoothers_[level] = smoother;
        }

        void init_memory(const Layout &l) override {
            Super::init_memory(l);
            VariableBoundSolverInterface::init_memory(l);
        }

        /**
         * @brief      The solve function for multigrid method.
         *
         * @param[in]  rhs   The right hand side.
         * @param      x   The initial guess.
         *
         */
        bool apply(const Vector &rhs, Vector &x_fine) override {
            Scalar r_norm, /*r0_norm,*/ diff_norm;
            SizeType it = 0;
            bool converged = false;
            bool ok = true;
            UTOPIA_UNUSED(ok);
            SizeType L = this->n_levels();

            memory.init(L);
            SizeType l = L - 1;

            memory.x[l] = x_fine;
            memory.rhs[l] = rhs;

            r_norm = norm2(memory.rhs[l] - level(l).A() * memory.x[l]);
            // r0_norm = r_norm;

            Vector x_old = memory.x[l];

            std::string mg_header_message = "Multigrid: " + std::to_string(L) + " levels";
            this->init_solver(mg_header_message, {" it. ", "|| r_N ||", "||x_old - x_new||"});

            if (this->verbose()) PrintInfo::print_iter_status(it, {r_norm, 1});
            it++;

            if (this->cycle_type() == FULL_CYCLE) {
                this->max_it(1);
            }

            while (!converged) {
                if (this->cycle_type() == MULTIPLICATIVE_CYCLE) {
                    ok = standard_cycle(l);
                    assert(ok);
                } else {
                    std::cout << "ERROR::UTOPIA_MG<< unknown MG type... \n";
                }

#ifdef CHECK_NUM_PRECISION_mode
                if (has_nan_or_inf(x_fine)) {
                    x_fine.set(0.0);
                    return true;
                }
#else
                // assert(!has_nan_or_inf(x));
#endif
                diff_norm = norm2(x_old - memory.x[l]);
                x_old = memory.x[l];

                r_norm = norm2(memory.rhs[l] - level(l).A() * memory.x[l]);

                // print iteration status on every iteration
                if (this->verbose()) PrintInfo::print_iter_status(it, {r_norm, diff_norm});

                // check convergence and print interation info
                converged = this->check_convergence(it, r_norm, diff_norm, 1);
                it++;
            }

            x_fine = memory.x[l];
            this->print_statistics(it);

            return true;
        }

        inline Level &level(const SizeType &l) { return this->levels_[l]; }

        /*=======================================================================================================================================
         =
         =========================================================================================================================================*/
    private:
        bool standard_cycle(const SizeType &l) {
            assert(memory.valid(this->n_levels()) && l < this->n_levels());

            ////////////////////////////////////
            if (l == 0) {
                coarse_solve(memory.rhs[l], memory.x[l]);
                return true;
            }
            ////////////////////////////////////

            // recursive call into mg
            for (SizeType k = 0; k < this->mg_type(); k++) {
                // presmoothing

                if (l == this->n_levels() - 1) {
                    // do fine level smoothing
                    smoothing_fine(l, memory.rhs[l], memory.x[l], this->pre_smoothing_steps(), true);
                } else if (l > 0) {
                    smoothing(l, memory.rhs[l], memory.x[l], this->pre_smoothing_steps());
                }

                memory.r[l] = memory.rhs[l] - level(l).A() * memory.x[l];
                // residual transfer
                this->transfer(l - 1).restrict(memory.r[l], memory.rhs[l - 1]);

                // BC conditions treatment ...
                if (this->must_generate_masks()) {
                    this->apply_mask(l - 1, memory.rhs[l - 1]);
                }

                if (empty(memory.x[l - 1])) {
                    memory.x[l - 1] = 0.0 * memory.rhs[l - 1];
                } else {
                    memory.x[l - 1].set(0.0);
                }

                standard_cycle(l - 1);

                // correction transfer
                this->transfer(l - 1).interpolate(memory.x[l - 1], memory.r[l]);

                memory.x[l] += memory.r[l];

                // postsmoothing
                if (l == this->n_levels() - 1) {
                    smoothing_fine(l, memory.rhs[l], memory.x[l], this->post_smoothing_steps(), false);
                } else if (l > 0) {
                    smoothing(l, memory.rhs[l], memory.x[l], this->post_smoothing_steps());
                }
            }

            return true;
        }

        /**
         * @brief      The function invokes smoothing.
         *
         * @param[in]  level The level we are at.
         * @param[in]  rhs   The right hand side.
         * @param      x     The current iterate.
         * @param[in]  nu    The number of smoothing steps.
         *
         * @return
         */
        inline bool smoothing(const SizeType l, const Vector &rhs, Vector &x, const SizeType &nu = 1) {
            smoothers_[l]->sweeps(nu);
            smoothers_[l]->smooth(rhs, x);
            return true;
        }

        inline bool smoothing_fine(const SizeType l,
                                   const Vector &rhs,
                                   Vector &x,
                                   const SizeType &nu = 1,
                                   const bool &pre_sm = false) {
            assert(l != this->n_levels() - 2);

            if (auto *GS_smoother = dynamic_cast<ProjectedGaussSeidelQR<Matrix, Vector> *>(smoothers_[l].get())) {
                GS_smoother->set_R(R_);
                GS_smoother->sweeps(nu);
                GS_smoother->set_box_constraints(this->get_box_constraints());
                GS_smoother->smooth(rhs, x);
            } else {
                utopia_error(
                    "MG_QR: requires ProjectedGaussSeidelQR to be the fine level "
                    "smoother ");
            }

            if (pre_sm) {
                if (auto *trunc_transfer =
                        dynamic_cast<IPRTruncatedTransfer<Matrix, Vector> *>(this->transfers_[l - 1].get())) {
                    auto *GS_smoother = dynamic_cast<ProjectedGaussSeidelQR<Matrix, Vector> *>(smoothers_[l].get());
                    const Vector &active_set = GS_smoother->get_active_set();

                    trunc_transfer->truncate_interpolation(active_set);
                    this->galerkin_assembly(this->get_operator());

                    for (std::size_t l = 1; l != smoothers_.size() - 1; ++l) {
                        smoothers_[l]->update(level(l).A_ptr());
                    }
                    coarse_solver_->update(level(0).A_ptr());
                }

                else {
                    utopia_error("MG_QR: requires IPRTruncatedTransfer for the finest level.");
                }
            }
            return true;
        }
        /**
         * @brief      The functions invokes coarse solve.
         *
         * @param[in]  rhs  The right hand side.
         * @param      x    The current iterate.
         *
         * @return
         */
        bool coarse_solve(const Vector &rhs, Vector &x) {
            if (!coarse_solver_->apply(rhs, x)) return false;
            assert(approxeq(level(0).A() * x, rhs, 1e-6));
            return true;
        }

        MultigridQR *clone() const override {
            return new MultigridQR(std::shared_ptr<Smoother>(smoothers_[smoothers_.size() - 1]->clone()),
                                   std::shared_ptr<Smoother>(smoothers_[1]->clone()),
                                   std::shared_ptr<Solver>(coarse_solver_->clone()),
                                   this->n_levels());
        }

    public:
        /**
         * @brief      Function changes direct solver needed for coarse grid solve.
         *
         * @param[in]  linear_solver  The linear solver.
         *
         * @return
         */
        bool change_coarse_solver(const std::shared_ptr<Solver> &coarse_solver) {
            coarse_solver_ = coarse_solver;
            return true;
        }

        void use_line_search(const bool val) { use_line_search_ = val; }

        bool use_line_search() const { return use_line_search_; }

        void set_R(const Matrix &R) { R_ = R; }

        void set_Q(const Matrix &Q) { Q_ = Q; }

        void set_QR(const Matrix &Q, const Matrix &R) {
            Q_ = Q;
            R_ = R;
        }

    protected:
        std::vector<SmootherPtr> smoothers_;
        std::shared_ptr<Solver> coarse_solver_;

    private:
        bool use_line_search_;
        Matrix Q_, R_;
    };

}  // namespace utopia

#endif  // UTOPIA_MULTIGRID_QR_HPP
