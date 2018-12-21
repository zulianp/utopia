#ifndef UTOPIA_PETSC_KSP_SOLVERS_HPP
#define UTOPIA_PETSC_KSP_SOLVERS_HPP

#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_petsc_ConjugateGradient.hpp"
#include "utopia_BiCGStab.hpp"

namespace utopia {
    //FIXME use superclass IterativeSolver instead of KSPSolver and compose with it
    template<typename Matrix, typename Vector>
    class BiCGStab<Matrix, Vector, PETSC> final : public KSPSolver<Matrix, Vector, PETSC> {
    public:
        using super = utopia::KSPSolver<Matrix, Vector, PETSC>;

        BiCGStab(const std::string &preconditioner = "jacobi")
        : KSPSolver<Matrix, Vector, PETSC>()
        {
            this->pc_type(preconditioner);
            this->ksp_type("bcgs");
        }

        void read(Input &is) override {
            super::read(is);
            super::ksp_type("bcgs");
        }

        virtual BiCGStab * clone() const override
        {
            return new BiCGStab(*this);
        }

    };

    //////////////////////////////////////////////////////////////////////////////////////////////////
    template<typename Matrix, typename Vector>
    class MINRES<Matrix, Vector, PETSC> final : public KSPSolver<Matrix, Vector, PETSC> {
    public:
        using super = utopia::KSPSolver<Matrix, Vector, PETSC>;

        MINRES(const std::string &preconditioner = "jacobi")
        : KSPSolver<Matrix, Vector, PETSC>()
        {
            this->pc_type(preconditioner);
            this->ksp_type("minres");
        }

        void read(Input &is) override {
            super::read(is);
            super::ksp_type("minres");
        }

        virtual MINRES * clone() const override
        {
            return new MINRES(*this);
        }
    };


    //////////////////////////////////////////////////////////////////////////////////////////////////
    template<typename Matrix, typename Vector>
    class SOR<Matrix, Vector, PETSC> final : public KSPSolver<Matrix, Vector, PETSC> {
    public:
        using super = utopia::KSPSolver<Matrix, Vector, PETSC>;
        
        SOR()
        : KSPSolver<Matrix, Vector, PETSC>()
        {
            this->pc_type("sor");
            this->ksp_type("richardson");
        }

        void read(Input &is) override {
            super::read(is);
            super::ksp_type("richardson");
            super::pc_type("sor");
        }

        virtual SOR * clone() const override
        {
            return new SOR(*this);
        }

    };
}
//////////////////////////////////////////////////////////////////////////////////////////////////

#endif //UTOPIA_PETSC_KSP_SOLVERS_HPP
