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
    : KSPSolver<Matrix, Vector, PETSC>(params), preconditioner_(preconditioner)
    {
      this->pc_type(preconditioner);
      this->ksp_type("bcgs");
    }

    void get_parameters(Parameters &params) const override {
        IterativeSolver<Matrix, Vector>::get_parameters(params);
        params.lin_solver_type("bcgs");
        params.preconditioner_type(preconditioner_.c_str());
    }

    void set_parameters(const Parameters params) override {
      Parameters params_copy = params;
      params_copy.lin_solver_type("bcgs");
      params_copy.preconditioner_type(preconditioner_.c_str());
      IterativeSolver<Matrix, Vector>::set_parameters(params_copy);
    }

    inline void pc_type(const std::string & preconditioner) override
    {
      preconditioner_ = preconditioner;
      KSPSolver<Matrix, Vector, PETSC>::pc_type(preconditioner_);
    }

    virtual BiCGStab * clone() const override 
    {
        return new BiCGStab(*this);
    }

  private:
    std::string preconditioner_;
  };

//////////////////////////////////////////////////////////////////////////////////////////////////
  template<typename Matrix, typename Vector> 
  class MINRES<Matrix, Vector, PETSC> : public KSPSolver<Matrix, Vector, PETSC> {
  public:
    MINRES(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
    : KSPSolver<Matrix, Vector, PETSC>(params), preconditioner_(preconditioner)
    {
      this->pc_type(preconditioner);
      this->ksp_type("minres");
    }

    void set_parameters(const Parameters params) override {
      Parameters params_copy = params;
      params_copy.lin_solver_type("minres");
      params_copy.preconditioner_type(preconditioner_.c_str());
      KSPSolver<Matrix, Vector, PETSC>::set_parameters(params_copy);
    }

    inline void pc_type(const std::string & preconditioner) override
    {
      preconditioner_ = preconditioner;
      KSPSolver<Matrix, Vector, PETSC>::pc_type(preconditioner_);
    }

    virtual MINRES * clone() const override 
    {
        return new MINRES(*this);
    }

  private:
    std::string preconditioner_;
  };

}
//////////////////////////////////////////////////////////////////////////////////////////////////

#endif //UTOPIA_PETSC_KSP_SOLVERS_HPP
