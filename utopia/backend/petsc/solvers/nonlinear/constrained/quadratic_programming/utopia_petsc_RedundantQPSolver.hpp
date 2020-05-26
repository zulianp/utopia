#ifndef UTOPIA_PETSC_REDUNDANT_QP_SOLVER_HPP
#define UTOPIA_PETSC_REDUNDANT_QP_SOLVER_HPP

#include "utopia_QPSolver.hpp"
#include "utopia_RedundantQPSolver.hpp"
#include "utopia_petsc_Redundant.hpp"

namespace utopia {

    template <>
    class RedundantQPSolver<PetscMatrix, PetscVector, PETSC> final
        : public OperatorBasedQPSolver<PetscMatrix, PetscVector> {
        using Traits = utopia::Traits<PetscVector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using IndexSet = typename Traits::IndexSet;
        using Comm = typename Traits::Communicator;
        using Layout = typename Traits::Layout;
        using MatrixLayout = typename Traits::MatrixLayout;

        using OperatorBasedQPSolver = utopia::OperatorBasedQPSolver<PetscMatrix, PetscVector>;
        using Super = OperatorBasedQPSolver;

    public:
        using Super::solve;
        using Super::update;

        RedundantQPSolver(const std::shared_ptr<OperatorBasedQPSolver> &qp_solver, const int n_sub_comm = 2)
            : qp_solver_(qp_solver) {
            red_.n_sub_comm(n_sub_comm);
        }

        void read(Input &in) override;

        ~RedundantQPSolver() override = default;
        RedundantQPSolver(const RedundantQPSolver &other);

        RedundantQPSolver *clone() const override { return new RedundantQPSolver(*this); }

        void update(const Operator<PetscVector> &A) override;
        bool solve(const Operator<PetscVector> &A, const PetscVector &rhs, PetscVector &sol) override;
        bool valid() const;

    private:
        std::shared_ptr<OperatorBasedQPSolver> qp_solver_;
        Redundant<PetscMatrix, PetscVector> red_;

        PetscMatrix redundant_matrix_;
        PetscVector redundant_sol_, redundant_rhs_;

        std::shared_ptr<PetscVector> redundant_lower_bound_;
        std::shared_ptr<PetscVector> redundant_upper_bound_;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_REDUNDANT_QP_SOLVER_HPP