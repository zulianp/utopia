#ifndef UTOPIA_ALGEBRAIC_MULTIGRID_HPP
#define UTOPIA_ALGEBRAIC_MULTIGRID_HPP

#include "utopia_Clonable.hpp"
#include "utopia_IPRTransfer.hpp"
#include "utopia_Multigrid.hpp"

namespace utopia {

    template <class Matrix>
    class MatrixAgglomerator : public Clonable, public Configurable {
    public:
        using Transfer = utopia::Transfer<Matrix, typename Traits<Matrix>::Vector>;

        virtual ~MatrixAgglomerator() = default;
        virtual std::shared_ptr<Transfer> create_transfer(const Matrix &in) = 0;
        virtual void read(Input &) override {}
        MatrixAgglomerator *clone() const override = 0;
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

        using Super::verbose;

        AlgebraicMultigrid *clone() const override { return new AlgebraicMultigrid(*this); }

        AlgebraicMultigrid(const AlgebraicMultigrid &other)
            : algo_(other.algo_),
              agglomerator_(std::shared_ptr<MatrixAgglomerator<Matrix>>(other.agglomerator_->clone())),
              n_levels_(other.n_levels_) {}

        void read(Input &in) override {
            Super::read(in);
            in.get("n_levels", n_levels_);
            algo_.read(in);

            if (agglomerator_) {
                in.get("agglomerator", *agglomerator_);
            }
        }

        bool apply(const Vector &rhs, Vector &sol) override {
            UTOPIA_TRACE_REGION_BEGIN("AlgebraicMultigrid::apply");
            bool ok = algo_.apply(rhs, sol);
            UTOPIA_TRACE_REGION_END("AlgebraicMultigrid::apply");
            return ok;
        }

        void update(const std::shared_ptr<const Matrix> &op) override {
            UTOPIA_TRACE_REGION_BEGIN("AlgebraicMultigrid::update");

            Super::update(op);

            //////////////////////////////////////////////////////////////////
            std::vector<std::shared_ptr<Transfer>> transfers(n_levels_ - 1);
            std::vector<std::shared_ptr<const Matrix>> matrices(n_levels_);
            matrices[n_levels_ - 1] = op;

            auto last_mat = op;

            for (SizeType l = n_levels_ - 2; l >= 0; --l) {
                auto A = std::make_shared<Matrix>();

                auto t = agglomerator_->create_transfer(*last_mat);

                auto temp_mat = std::make_shared<Matrix>();
                t->restrict(*last_mat, *temp_mat);

                transfers[l] = t;
                matrices[l] = temp_mat;
                last_mat = temp_mat;
            }

            //////////////////////////////////////////////////////////////////
            algo_.set_transfer_operators(transfers);
            algo_.set_linear_operators(matrices);
            algo_.set_perform_galerkin_assembly(false);
            algo_.update();

            UTOPIA_TRACE_REGION_END("AlgebraicMultigrid::update");
        }

        void verbose(const bool &val) override {
            Super::verbose(val);
            algo_.verbose(val);
        }

        AlgebraicMultigrid(const std::shared_ptr<Smoother> &smoother,
                           const std::shared_ptr<Solver> &coarse_solver,
                           const std::shared_ptr<MatrixAgglomerator<Matrix>> &agglomerator)
            : algo_(smoother, coarse_solver), agglomerator_(agglomerator) {}

    private:
        Multigrid<Matrix, Vector> algo_;
        std::shared_ptr<MatrixAgglomerator<Matrix>> agglomerator_;
        SizeType n_levels_{2};
    };

}  // namespace utopia

#endif  // UTOPIA_ALGEBRAIC_MULTIGRID_HPP
