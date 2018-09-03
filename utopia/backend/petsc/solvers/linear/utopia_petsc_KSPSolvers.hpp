#ifndef UTOPIA_PETSC_KSP_SOLVERS_HPP
#define UTOPIA_PETSC_KSP_SOLVERS_HPP

#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_petsc_ConjugateGradient.hpp"

namespace utopia {
    //FIXME use superclass IterativeSolver instead of KSPSolver and compose with it
    template<typename Matrix, typename Vector>
    class BiCGStab<Matrix, Vector, PETSC> : public KSPSolver<Matrix, Vector, PETSC> {
    public:
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

        virtual BiCGStab * clone() const override
        {
            return new BiCGStab(*this);
        }

    };

    //////////////////////////////////////////////////////////////////////////////////////////////////
    template<typename Matrix, typename Vector>
    class MINRES<Matrix, Vector, PETSC> : public KSPSolver<Matrix, Vector, PETSC> {
    public:
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

        virtual MINRES * clone() const override
        {
            return new MINRES(*this);
        }
    };


    //////////////////////////////////////////////////////////////////////////////////////////////////
    template<typename Matrix, typename Vector>
    class SOR<Matrix, Vector, PETSC> : public KSPSolver<Matrix, Vector, PETSC> {
    public:
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

        virtual SOR * clone() const override
        {
            return new SOR(*this);
        }

    };
}
//////////////////////////////////////////////////////////////////////////////////////////////////

#endif //UTOPIA_PETSC_KSP_SOLVERS_HPP
