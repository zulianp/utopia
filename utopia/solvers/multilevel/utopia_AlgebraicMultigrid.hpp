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
        virtual std::shared_ptr<Transfer> create_truncated_transfer(const Matrix &in) { return create_transfer(in); }
        virtual void read(Input &) override {}
        MatrixAgglomerator *clone() const override = 0;
    };

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

        using Super::max_it;
        using Super::verbose;

        AlgebraicMultigrid *clone() const override { return new AlgebraicMultigrid(*this); }

        AlgebraicMultigrid(const AlgebraicMultigrid &other)
            : algorithm_(other.algorithm_),
              agglomerator_(std::shared_ptr<MatrixAgglomerator<Matrix>>(other.agglomerator_->clone())),
              n_levels_(other.n_levels_) {}

        void read(Input &in) override {
            Super::read(in);
            in.get("n_levels", n_levels_);
            algorithm_.read(in);

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

        void update(const std::shared_ptr<const Matrix> &op) override {
            UTOPIA_TRACE_REGION_BEGIN("AlgebraicMultigrid::update");

            Super::update(op);

            AlgebraicMultigridBuilder<Matrix>::build(n_levels_, op, *agglomerator_, algorithm_);
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

        AlgebraicMultigrid(const std::shared_ptr<Smoother> &smoother,
                           const std::shared_ptr<Solver> &coarse_solver,
                           const std::shared_ptr<MatrixAgglomerator<Matrix>> &agglomerator)
            : algorithm_(smoother, coarse_solver), agglomerator_(agglomerator) {}

        inline const Multigrid<Matrix, Vector> &algorithm() const { return algorithm_; }
        inline Multigrid<Matrix, Vector> &algorithm() { return algorithm_; }

    private:
        Multigrid<Matrix, Vector> algorithm_;
        std::shared_ptr<MatrixAgglomerator<Matrix>> agglomerator_;
        SizeType n_levels_{2};
    };

}  // namespace utopia

#endif  // UTOPIA_ALGEBRAIC_MULTIGRID_HPP
