#ifndef UTOPIA_PETSC_SCHUR_COMPLEMENT_HPP
#define UTOPIA_PETSC_SCHUR_COMPLEMENT_HPP

#include "utopia_Operator.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"

#include "utopia_Input.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_petsc_Traits.hpp"

#include <memory>

namespace utopia {

    template <class Matrix>
    class SchurComplement {};

    template <>
    class SchurComplement<PetscMatrix> : public Operator<PetscVector>, public Configurable {
    public:
        using IndexSet = typename Traits<PetscMatrix>::IndexSet;
        using IndexArray = typename Traits<PetscMatrix>::IndexArray;
        using Communicator = typename Traits<PetscMatrix>::Communicator;
        using SizeType = typename Traits<PetscMatrix>::SizeType;
        using Scalar = typename Traits<PetscMatrix>::Scalar;

        bool apply(const PetscVector &x, PetscVector &y) const override;

        bool apply_righthand_side(const PetscVector &rhs, PetscVector &out);
        bool initialize_from_selection(const PetscMatrix &matrix, const IndexArray &eliminated_dofs);
        bool finalize(const PetscVector &rhs, const PetscVector &x_restricted, PetscVector &x);

        void set_solver(const std::shared_ptr<LinearSolver<PetscMatrix, PetscVector>> &solver);

        void read(Input &in) override;

        Size size() const override;
        Size local_size() const override;
        Communicator &comm() override;
        const Communicator &comm() const override;

        ~SchurComplement();
        SchurComplement();

        class Impl;
        std::unique_ptr<Impl> impl_;
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_SCHUR_COMPLEMENT_HPP
