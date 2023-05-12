#ifndef UTOPIA_PETSC_REBALANCED_SOLVER_HPP
#define UTOPIA_PETSC_REBALANCED_SOLVER_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_PARMETIS

#include "utopia_LinearSolver.hpp"

#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

namespace utopia {
    //
    class RebalancedSolver : public LinearSolver<PetscMatrix, PetscVector> {
    public:
        using Super = utopia::LinearSolver<PetscMatrix, PetscVector>;

        RebalancedSolver *clone() const override;
        bool apply(const PetscVector &rhs, PetscVector &sol) override;
        void update(const std::shared_ptr<const PetscMatrix> &op) override;
        void read(Input &in) override;

        RebalancedSolver();
        ~RebalancedSolver();

        void set_solver(const std::shared_ptr<LinearSolver<PetscMatrix, PetscVector>> &solver);

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };
}  // namespace utopia

#endif

#endif  // UTOPIA_PETSC_REBALANCED_SOLVER_HPP
