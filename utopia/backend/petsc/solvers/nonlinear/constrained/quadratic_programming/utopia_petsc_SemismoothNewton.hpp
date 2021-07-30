#ifndef UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP
#define UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP

#include <memory>
#include "utopia_QPSolver.hpp"
#include "utopia_SemismoothNewton.hpp"
#include "utopia_petsc.hpp"
#include "utopia_petsc_impl.hpp"

#include "utopia_petsc_Factorization.hpp"

namespace utopia {

    template <>
    class SemismoothNewton<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL> final
        : public QPSolver<PetscMatrix, PetscVector> {
        using Scalar = utopia::Traits<PetscVector>::Scalar;
        using SizeType = utopia::Traits<PetscVector>::SizeType;
        typedef utopia::LinearSolver<PetscMatrix, PetscVector> LinearSolver;
        typedef utopia::QPSolver<PetscMatrix, PetscVector> Super;
        using BoxConstraints = utopia::BoxConstraints<PetscVector>;

    public:
        SemismoothNewton(const std::shared_ptr<LinearSolver> &linear_solver =
                             std::make_shared<Factorization<PetscMatrix, PetscVector>>());
        ~SemismoothNewton() override;
        SemismoothNewton *clone() const override;
        bool apply(const PetscVector &rhs, PetscVector &sol) override;
        void read(Input &in) override;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP
