#include "utopia_petsc_SemismoothNewton.hpp"
#include "utopia_QuadraticFunction.hpp"
#include "utopia_petsc_SNES.hpp"

#include "utopia_make_unique.hpp"
#include <cassert>

namespace utopia {

    template<class Matrix, class Vector>
    class SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL>::Impl {
    public:
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;

        Impl(const std::shared_ptr<LinearSolver> &linear_solver)
        : snes(utopia::make_unique<SNESSolver<Matrix, Vector>>(linear_solver))
        {
            snes->set_snes_type(SNESVINEWTONSSLS);
            snes->set_line_search_type(SNESLINESEARCHBASIC);
        }

        std::unique_ptr<SNESSolver<Matrix, Vector>> snes;
    };


    template<class Matrix, class Vector>
    SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL>::SemismoothNewton(const std::shared_ptr<LinearSolver> &linear_solver)
    : impl_(utopia::make_unique<Impl>(linear_solver))
    { }

    template <class Matrix, class Vector>
    SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL>::~SemismoothNewton() = default;

    template<class Matrix, class Vector>
    SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL> * SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL>::clone() const {
        //exception-safe cloning
        auto cloned_snes    = std::unique_ptr<SNESSolver<Matrix, Vector>>(impl_->snes->clone());
        auto cloned         = utopia::make_unique<SemismoothNewton>();
        cloned->impl_->snes = std::move(cloned_snes);
        return cloned.release();
    }

    template<class Matrix, class Vector>
    bool SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL>::apply(const Vector &rhs, Vector &sol)
    {
        assert(this->has_operator());

        impl_->snes->atol(this->atol());
        impl_->snes->stol(this->stol());
        impl_->snes->rtol(this->rtol());
        impl_->snes->max_it(this->max_it());

        impl_->snes->set_box_constraints(this->get_box_constraints());
        //(JUST TO BE SURE) FIXME find out why this changes
        impl_->snes->set_snes_type(SNESVINEWTONSSLS);

        QuadraticFunction<Matrix, Vector> fun(
            std::make_shared<Matrix>(*this->get_operator()),
            std::make_shared<Vector>(rhs)
        );

        return impl_->snes->solve(fun, sol);
    }

    template<class Matrix, class Vector>
    void SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL>::read(Input &in)
    {
        impl_->snes->read(in);
    }

}

