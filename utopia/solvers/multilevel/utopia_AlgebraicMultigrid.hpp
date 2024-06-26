#ifndef UTOPIA_ALGEBRAIC_MULTIGRID_HPP
#define UTOPIA_ALGEBRAIC_MULTIGRID_HPP

#include "utopia_Clonable.hpp"
#include "utopia_IPRTransfer.hpp"
#include "utopia_Multigrid.hpp"

#include "utopia_Agglomerate.hpp"
#include "utopia_BlockAgglomerate.hpp"
#include "utopia_ILU.hpp"
#include "utopia_MatrixAgglomerator.hpp"

namespace utopia {

    template <class Matrix>
    class AlgebraicMultigridBuilder {
    public:
        using Vector = typename Traits<Matrix>::Vector;
        using Transfer = utopia::Transfer<Matrix, Vector>;

        static void build(const int n_levels,
                          const std::shared_ptr<const Matrix> &op,
                          MatrixAgglomerator<Matrix> &agglomerator,
                          Multigrid<Matrix, Vector> &mg) {
            UTOPIA_TRACE_REGION_BEGIN("AlgebraicMultigridBuilder::build");

            //////////////////////////////////////////////////////////////////
            std::vector<std::shared_ptr<Transfer>> transfers(n_levels - 1);
            std::vector<std::shared_ptr<const Matrix>> matrices(n_levels);
            matrices[n_levels - 1] = op;

            auto last_mat = op;

            for (SizeType l = n_levels - 2; l >= 0; --l) {
                auto A = std::make_shared<Matrix>();

                auto t = agglomerator.create_transfer(*last_mat);

                auto temp_mat = std::make_shared<Matrix>();
                t->restrict(*last_mat, *temp_mat);

                transfers[l] = t;
                matrices[l] = temp_mat;
                last_mat = temp_mat;
            }

            //////////////////////////////////////////////////////////////////
            mg.set_transfer_operators(transfers);
            mg.set_linear_operators(matrices);
            mg.set_perform_galerkin_assembly(false);
            mg.update();

            UTOPIA_TRACE_REGION_END("AlgebraicMultigridBuilder::build");
        }
    };

    template <class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class AlgebraicMultigrid final : public IterativeSolver<Matrix, Vector> {
    public:
        using Super = utopia::IterativeSolver<Matrix, Vector>;
        using Solver = utopia::LinearSolver<Matrix, Vector>;
        using IterativeSolver = utopia::IterativeSolver<Matrix, Vector>;
        using Smoother = utopia::IterativeSolver<Matrix, Vector>;
        using Level = utopia::Level<Matrix, Vector>;
        using Transfer = utopia::Transfer<Matrix, Vector>;
        using SmootherPtr = std::shared_ptr<Smoother>;
        using Scalar = typename Traits<Matrix>::Scalar;
        using SizeType = typename Traits<Matrix>::SizeType;

        using ILU = utopia::ILU<Matrix, Vector>;

        using Super::max_it;
        using Super::verbose;

        AlgebraicMultigrid *clone() const override { return new AlgebraicMultigrid(*this); }

        AlgebraicMultigrid(const AlgebraicMultigrid &other)
            : IterativeSolver(other),
              algorithm_(other.algorithm_),
              agglomerator_(std::shared_ptr<MatrixAgglomerator<Matrix>>(other.agglomerator_->clone())),
              n_levels_(other.n_levels_) {}

        inline void set_n_levels(const int n_levels) { n_levels_ = n_levels; }

        inline void set_block_size(const int block_size) { block_size_ = block_size; }

        void ensure_agglomerator() {
            if (!agglomerator_) {
#ifdef UTOPIA_ENABLE_PETSC
                if (block_size_ == 2) {
                    agglomerator_ = std::make_shared<BlockAgglomerate<Matrix, 2>>();
                } else if (block_size_ == 3) {
                    agglomerator_ = std::make_shared<BlockAgglomerate<Matrix, 3>>();
                } else if (block_size_ == 4) {
                    agglomerator_ = std::make_shared<BlockAgglomerate<Matrix, 4>>();
                } else
#endif  // UTOPIA_ENABLE_PETSC

                {
                    agglomerator_ = std::make_shared<Agglomerate<Matrix>>();
                }
            }
        }

        void read(Input &in) override {
            Super::read(in);
            in.get("n_levels", n_levels_);
            in.get("block_size", block_size_);

            algorithm_.read(in);

            ensure_agglomerator();

            if (agglomerator_) {
                in.get("agglomerator", *agglomerator_);
            }
        }

        bool apply(const Vector &rhs, Vector &sol) override {
            UTOPIA_TRACE_REGION_BEGIN("AlgebraicMultigrid::apply");
            bool ok = algorithm_.apply(rhs, sol);
            UTOPIA_TRACE_REGION_END("AlgebraicMultigrid::apply");
            return ok;
        }

        bool smooth(const Vector &rhs, Vector &x) override {
            SizeType temp = this->max_it();
            this->max_it(this->sweeps());
            this->apply(rhs, x);
            this->max_it(temp);
            return true;
        }

        void update(const std::shared_ptr<const Matrix> &op) override {
            UTOPIA_TRACE_REGION_BEGIN("AlgebraicMultigrid::update");

            Super::update(op);

            ensure_agglomerator();
            AlgebraicMultigridBuilder<Matrix>::build(n_levels_, op, *agglomerator_, algorithm_);

            algorithm_.adjust_memory();
            UTOPIA_TRACE_REGION_END("AlgebraicMultigrid::update");
        }

        void verbose(const bool &val) override {
            Super::verbose(val);
            algorithm_.verbose(val);
        }

        void max_it(const SizeType &max_it_in) override {
            Super::max_it(max_it_in);
            algorithm_.max_it(max_it_in);
        }

        void set_agglomerator(const std::shared_ptr<MatrixAgglomerator<Matrix>> &agglomerator) {
            agglomerator_ = agglomerator;
        }

        AlgebraicMultigrid(const std::shared_ptr<Smoother> &smoother = std::make_shared<ILU>(),
                           const std::shared_ptr<Solver> &coarse_solver = std::make_shared<ILU>(),
                           const std::shared_ptr<MatrixAgglomerator<Matrix>> &agglomerator = nullptr)
            : algorithm_(smoother, coarse_solver), agglomerator_(agglomerator) {}

        inline const Multigrid<Matrix, Vector> &algorithm() const { return algorithm_; }
        inline Multigrid<Matrix, Vector> &algorithm() { return algorithm_; }

    private:
        Multigrid<Matrix, Vector> algorithm_;
        std::shared_ptr<MatrixAgglomerator<Matrix>> agglomerator_;
        SizeType n_levels_{2};
        int block_size_{1};
    };  // namespace utopia

}  // namespace utopia

#endif  // UTOPIA_ALGEBRAIC_MULTIGRID_HPP
