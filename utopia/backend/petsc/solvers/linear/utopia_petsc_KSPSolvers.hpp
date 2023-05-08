#ifndef UTOPIA_PETSC_KSP_SOLVERS_HPP
#define UTOPIA_PETSC_KSP_SOLVERS_HPP

#include "utopia_BiCGStab.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_ConjugateGradient.hpp"
#include "utopia_petsc_KSPSolver.hpp"

namespace utopia {
    // FIXME use superclass IterativeSolver instead of KSPSolver and compose with it
    template <typename Matrix, typename Vector>
    class BiCGStab<Matrix, Vector, PETSC> final : public KSPSolver<Matrix, Vector, PETSC> {
    public:
        using super = utopia::KSPSolver<Matrix, Vector, PETSC>;

        BiCGStab(const std::string &preconditioner = PCJACOBI) : KSPSolver<Matrix, Vector, PETSC>() {
            this->pc_type(preconditioner);
            this->ksp_type(KSPBCGS);
        }

        void read(Input &is) override {
            super::read(is);
            super::ksp_type(KSPBCGS);
        }

        BiCGStab *clone() const override { return new BiCGStab(*this); }
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////
    template <typename Matrix, typename Vector>
    class MINRES<Matrix, Vector, PETSC> final : public KSPSolver<Matrix, Vector, PETSC> {
    public:
        using super = utopia::KSPSolver<Matrix, Vector, PETSC>;

        MINRES(const std::string &preconditioner = PCJACOBI) : KSPSolver<Matrix, Vector, PETSC>() {
            this->pc_type(preconditioner);
            this->ksp_type(KSPMINRES);
        }

        void read(Input &is) override {
            super::read(is);
            super::ksp_type(KSPMINRES);
        }

        MINRES *clone() const override { return new MINRES(*this); }
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////
    template <typename Matrix, typename Vector>
    class SOR<Matrix, Vector, PETSC> final : public KSPSolver<Matrix, Vector, PETSC> {
    public:
        using super = utopia::KSPSolver<Matrix, Vector, PETSC>;

        SOR() : KSPSolver<Matrix, Vector, PETSC>() {
            this->pc_type(PCSOR);
            this->ksp_type(KSPRICHARDSON);
        }

        void read(Input &is) override {
            super::read(is);
            super::ksp_type(KSPRICHARDSON);
            super::pc_type(PCSOR);
        }

        SOR *clone() const override { return new SOR(*this); }
    };
}  // namespace utopia
//////////////////////////////////////////////////////////////////////////////////////////////////

#endif  // UTOPIA_PETSC_KSP_SOLVERS_HPP
