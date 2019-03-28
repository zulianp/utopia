#include "utopia_polymorphic_LinearSolver.hpp"
#include "utopia_LinearSolverFactory.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    void PolymorphicLinearSolver<Matrix, Vector>::read(Input &in)
    {
        Super::read(in);

        std::string backend = UTOPIA_BACKEND(Vector).info().get_name();
        std::string type  = "gmres";

        in.get("backend", backend);
        in.get("type",    type);

        impl_ = LinearSolverFactory<Matrix, Vector, PETSC>::new_linear_solver(type);
        impl_->read(in);
    }

    template<class Matrix, class Vector>
    PolymorphicLinearSolver<Matrix, Vector>::PolymorphicLinearSolver()
    : impl_(LinearSolverFactory<Matrix, Vector, PETSC>::default_linear_solver())
    { }

    template<class Matrix, class Vector>
    PolymorphicLinearSolver<Matrix, Vector>::~PolymorphicLinearSolver()
    { }

    template<class Matrix, class Vector>
    PolymorphicLinearSolver<Matrix, Vector> * PolymorphicLinearSolver<Matrix, Vector>::clone() const
    {
        auto cloned = utopia::make_unique<PolymorphicLinearSolver>();

        if(this->impl_) {
            cloned->impl_ = std::unique_ptr<LinearSolver<Matrix, Vector>>(this->impl_->clone());
        }

        return cloned.release();
    }

    template<class Matrix, class Vector>
    bool PolymorphicLinearSolver<Matrix, Vector>::apply(const Vector &rhs, Vector &sol)
    {
        assert(static_cast<bool>(impl_));
        if(!impl_) return false;

        return impl_->apply(rhs, sol);
    }


    template<class Matrix, class Vector>
    void PolymorphicLinearSolver<Matrix, Vector>::update(const std::shared_ptr<const Matrix> &mat)
    {
        assert(static_cast<bool>(impl_));

        Super::update(mat);
        impl_->update(mat);
    }

}

