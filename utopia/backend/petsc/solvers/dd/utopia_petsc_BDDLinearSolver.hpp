#ifndef UTOPIA_PETSC_BDD_LINEAR_SOLVER_HPP
#define UTOPIA_PETSC_BDD_LINEAR_SOLVER_HPP

#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_PreconditionedSolver.hpp"

#include <memory>

namespace utopia {

    template <typename Matrix, typename Vector>
    class BDDLinearSolver : public PreconditionedSolver<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Preconditioner = utopia::Preconditioner<Vector>;
        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;
        using PreconditionedSolver = utopia::PreconditionedSolver<Matrix, Vector>;
        using Super = PreconditionedSolver;
        using Layout = typename Super::Layout;
        using MatrixFreeLinearSolver = utopia::MatrixFreeLinearSolver<Vector>;

        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "only works with petsc types");

        void read(Input &in) override;

        void init_memory(const Layout &layout) override;

        void print_usage(std::ostream &os) const override;

        bool apply(const Vector &b, Vector &x) override;

        void update(const std::shared_ptr<const Matrix> &A) override;

        void set_inner_solver(const std::shared_ptr<MatrixFreeLinearSolver> &solver);

        void verbose(const bool &val) override;
        void atol(const Scalar &val) override;
        void rtol(const Scalar &val) override;
        void stol(const Scalar &val) override;

        BDDLinearSolver *clone() const override;

        void copy(const BDDLinearSolver &other);
        BDDLinearSolver(const BDDLinearSolver &other);

        ~BDDLinearSolver();
        BDDLinearSolver();

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

    template <typename Matrix, typename Vector>
    class Traits<BDDLinearSolver<Matrix, Vector>> : public Traits<Matrix> {};
}  // namespace utopia

#endif  // UTOPIA_PETSC_BDD_LINEAR_SOLVER_HPP