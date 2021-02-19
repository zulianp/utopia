#ifndef UTOPIA_ALGEBRAIC_MULTIGRID_HPP
#define UTOPIA_ALGEBRAIC_MULTIGRID_HPP

#include "utopia_Clonable.hpp"
#include "utopia_IPRTransfer.hpp"
#include "utopia_Multigrid.hpp"

namespace utopia {

    template <class Matrix>
    class MatrixAgglomerator : public Clonable, public Configurable {
    public:
        virtual ~MatrixAgglomerator() = default;
        virtual void create_prolongator(const Matrix &in, Matrix &prolongator) = 0;
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

        bool apply(const Vector &rhs, Vector &sol) override { return algo_.apply(rhs, sol); }

        void update(const std::shared_ptr<const Matrix> &op) override {
            Super::update(op);

            //////////////////////////////////////////////////////////////////
            std::vector<std::shared_ptr<Transfer>> transfers(n_levels_ - 1);
            std::vector<std::shared_ptr<const Matrix>> matrices(n_levels_);
            matrices[n_levels_ - 1] = op;

            auto last_mat = op;

            for (SizeType l = n_levels_ - 2; l >= 0; --l) {
                auto A = std::make_shared<Matrix>();
                auto P = std::make_shared<Matrix>();
                agglomerator_->create_prolongator(*last_mat, *P);

                auto temp_mat = std::make_shared<Matrix>(transpose(*P) * (*last_mat) * (*P));

                // algo_.fix_semidefinite_operator(*temp_mat);

                transfers[l] = std::make_shared<IPRTransfer<Matrix, Vector>>(P);
                matrices[l] = temp_mat;
                last_mat = temp_mat;

                if (this->verbose()) {
                    auto coarsening_factor = P->rows() / Scalar(P->cols());
                    if (op->comm().rank() == 0) {
                        out() << "coarsening_factor(" << l << "): " << coarsening_factor << "\n";
                    }
                }
            }

            //////////////////////////////////////////////////////////////////
            algo_.set_transfer_operators(transfers);
            algo_.set_linear_operators(matrices);
            algo_.set_perform_galerkin_assembly(false);
            algo_.update();
        }

        void verbose(const bool &val) override {
            Super::verbose(val);
            algo_.verbose(false);
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
