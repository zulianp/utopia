#include "utopia_petsc_SemismoothNewton.hpp"

#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector_impl.hpp"

#include "utopia_QuadraticFunction.hpp"
#include "utopia_petsc_SNES.hpp"
#include "utopia_petsc_SemismoothNewton.hpp"

#include <cassert>
#include "utopia_make_unique.hpp"

namespace utopia {

    class SemismoothNewton<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>::Impl {
    public:
        typedef utopia::LinearSolver<PetscMatrix, PetscVector> LinearSolver;

        Impl(const std::shared_ptr<LinearSolver> &linear_solver)
            : snes(utopia::make_unique<SNESSolver<PetscMatrix, PetscVector>>(linear_solver)) {
            snes->set_snes_type(SNESVINEWTONSSLS);
            snes->set_line_search_type(SNESLINESEARCHBASIC);
        }

        std::unique_ptr<SNESSolver<PetscMatrix, PetscVector>> snes;
    };

    SemismoothNewton<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>::SemismoothNewton(
        const std::shared_ptr<LinearSolver> &linear_solver)
        : impl_(utopia::make_unique<Impl>(linear_solver)) {}

    SemismoothNewton<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>::~SemismoothNewton() = default;

    SemismoothNewton<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>
        *SemismoothNewton<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>::clone() const {
        // exception-safe cloning
        auto cloned_snes = std::unique_ptr<SNESSolver<PetscMatrix, PetscVector>>(impl_->snes->clone());
        auto cloned = utopia::make_unique<SemismoothNewton>();
        cloned->impl_->snes = std::move(cloned_snes);
        return cloned.release();
    }

    bool SemismoothNewton<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>::apply(const PetscVector &rhs,
                                                                               PetscVector &sol) {
        assert(this->has_operator());

        impl_->snes->atol(this->atol());
        impl_->snes->stol(this->stol());
        impl_->snes->rtol(this->rtol());
        impl_->snes->max_it(this->max_it());

        impl_->snes->set_box_constraints(this->get_box_constraints());
        //(JUST TO BE SURE) FIXME find out why this changes
        impl_->snes->set_snes_type(SNESVINEWTONSSLS);

        QuadraticFunction<PetscMatrix, PetscVector> fun(std::make_shared<PetscMatrix>(*this->get_operator()),
                                                        std::make_shared<PetscVector>(rhs));

        return impl_->snes->solve(fun, sol);
    }

    void SemismoothNewton<PetscMatrix, PetscVector, PETSC_EXPERIMENTAL>::read(Input &in) { impl_->snes->read(in); }

}  // namespace utopia
