#ifndef UTOPIA_PETSC_KSP_SOLVERS_HPP
#define UTOPIA_PETSC_KSP_SOLVERS_HPP

#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_petsc_ConjugateGradient.hpp"

namespace utopia {
    //FIXME use superclass IterativeSolver instead of KSPSolver and compose with it
    template<typename Matrix, typename Vector>
    class BiCGStab<Matrix, Vector, PETSC> final : public KSPSolver<Matrix, Vector, PETSC> {
    public:
        using super = utopia::KSPSolver<Matrix, Vector, PETSC>;

        BiCGStab(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
        : KSPSolver<Matrix, Vector, PETSC>(params)
        {
            this->pc_type(preconditioner);
            this->ksp_type("bcgs");
        }

        void get_parameters(Parameters &params) const override {
            IterativeSolver<Matrix, Vector>::get_parameters(params);
            params.lin_solver_type("bcgs");
            params.preconditioner_type(this->pc_type().c_str());
        }

        void set_parameters(const Parameters params) override {
            Parameters params_copy = params;
            params_copy.lin_solver_type("bcgs");
            params_copy.preconditioner_type(this->pc_type().c_str());
            IterativeSolver<Matrix, Vector>::set_parameters(params_copy);
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

        MINRES(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
        : KSPSolver<Matrix, Vector, PETSC>(params)
        {
            this->pc_type(preconditioner);
            this->ksp_type("minres");
        }

        void set_parameters(const Parameters params) override {
            Parameters params_copy = params;
            params_copy.lin_solver_type("minres");
            params_copy.preconditioner_type(this->pc_type().c_str());
            KSPSolver<Matrix, Vector, PETSC>::set_parameters(params_copy);
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
        
        SOR(const Parameters params = Parameters())
        : KSPSolver<Matrix, Vector, PETSC>(params)
        {
            this->pc_type("sor");
            this->ksp_type("richardson");
        }

        void set_parameters(const Parameters params) override {
            Parameters params_copy = params;
            params_copy.lin_solver_type("richardson");
            params_copy.preconditioner_type("sor");
            KSPSolver<Matrix, Vector, PETSC>::set_parameters(params_copy);
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
