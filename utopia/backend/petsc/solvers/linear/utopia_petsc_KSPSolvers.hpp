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

  private:
    std::string preconditioner_;
  };

//////////////////////////////////////////////////////////////////////////////////////////////////
  // template<typename Matrix, typename Vector> 
  // class ConjugateGradient<Matrix, Vector, PETSC> : public KSPSolver<Matrix, Vector, PETSC> {
  // public:
  //   ConjugateGradient(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
  //   : KSPSolver<Matrix, Vector, PETSC>(params), preconditioner_(preconditioner)
  //   {
  //     this->pc_type(preconditioner_);
  //     this->ksp_type("cg");
  //   }

  //   void set_parameters(const Parameters params) override {
  //     Parameters params_copy = params;
  //     params_copy.lin_solver_type("cg");
  //     params_copy.preconditioner_type(preconditioner_.c_str());
  //     KSPSolver<Matrix, Vector, PETSC>::set_parameters(params_copy);
  //   }

  //   inline void pc_type(const std::string & preconditioner) override
  //   {
  //     preconditioner_ = preconditioner;
  //     KSPSolver<Matrix, Vector, PETSC>::pc_type(preconditioner_);
  //   }

  // private:
  //   std::string preconditioner_;
  // };

//////////////////////////////////////////////////////////////////////////////////////////////////
  template<typename Matrix, typename Vector> 
  class GMRES<Matrix, Vector, PETSC> : public KSPSolver<Matrix, Vector, PETSC> {
  public:
    GMRES(const Parameters params = Parameters(), const std::string &preconditioner = "jacobi")
    : KSPSolver<Matrix, Vector, PETSC>(params), preconditioner_(preconditioner)
    {
      this->pc_type(preconditioner);
      this->ksp_type("gmres");
    }

    void set_parameters(const Parameters params) override {
      Parameters params_copy = params;
      params_copy.lin_solver_type("gmres");
      params_copy.preconditioner_type(preconditioner_.c_str());
      KSPSolver<Matrix, Vector, PETSC>::set_parameters(params_copy);
    }

    inline void pc_type(const std::string & preconditioner) override
    {
      preconditioner_ = preconditioner;
      KSPSolver<Matrix, Vector, PETSC>::pc_type(preconditioner_);
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


  private:
    std::string preconditioner_;
  };

//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// solver /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

  // template<typename Matrix, typename Vector> 
  // class Multigrid<Matrix, Vector, PETSC> :  public KSPSolver<Matrix, Vector, PETSC>,
  //                                           public MultiLevelBase<Matrix, Vector>  
  // {
  // public:
  //   Multigrid(const std::shared_ptr<Smoother<Matrix, Vector> > &smoother = std::shared_ptr<Smoother<Matrix, Vector> >(), 
  //             const std::shared_ptr<LinearSolver<Matrix, Vector> > &direct_solver = std::shared_ptr<LinearSolver<Matrix, Vector> >(),
  //             const Parameters params = Parameters(), 
  //             const std::string &preconditioner = "mg")
  //   : KSPSolver<Matrix, Vector, PETSC>(params), preconditioner_(preconditioner)
  //   {
  //     this->pc_type(preconditioner_);
  //     this->ksp_type("richardson");
  //     std::cout<<"UTOPIA:MG:PETSC.......... \n"; 
  //   }

  //   void set_parameters(const Parameters params) override {
  //     Parameters params_copy = params;
  //     params_copy.lin_solver_type("richardson");
  //     params_copy.preconditioner_type(preconditioner_.c_str());
  //     KSPSolver<Matrix, Vector, PETSC>::set_parameters(params_copy);
  //   }

  //   inline void pc_type(const std::string & preconditioner) override
  //   {
  //     preconditioner_ = preconditioner;
  //     KSPSolver<Matrix, Vector, PETSC>::pc_type(preconditioner_);
  //   }

  //   bool change_direct_solver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &linear_solver = std::shared_ptr<LinearSolver<Matrix, Vector>>())
  //   {
  //       std::cout<<"utopia::Multigrid::Petsc_backend::change_direct_solver --- not suported at the moment --- \n"; 
  //       return true; 
  //   }

  //   bool change_smoother(const std::shared_ptr<Smoother<Matrix, Vector> > &smoother = std::shared_ptr<Smoother<Matrix, Vector> >())
  //   {
  //       std::cout<<"utopia::Multigrid::Petsc_backend::change_smoother --- not suported at the moment --- \n"; 
  //       //PCMGGetSmoother(pc, l-1, &cksp); - one by one... 
  //       return true; 
  //   }

  //   void pre_smoothing_steps(const SizeType & pre_smoothing_steps_in ) 
  //   { 
  //       _pre_smoothing_steps = pre_smoothing_steps_in; 
  //       PC pc; 
  //       KSPGetPC(this->ksp, &pc);
  //       PCMGSetNumberSmoothDown(pc, pre_smoothing_steps_in);
  //   }; 

  //   void post_smoothing_steps(const SizeType & post_smoothing_steps_in )
  //   { 
  //       _post_smoothing_steps = post_smoothing_steps_in; 
  //       PC pc; 
  //       KSPGetPC(this->ksp, &pc);
  //       PCMGSetNumberSmoothUp(pc, post_smoothing_steps_in);
  //   }; 

  //   void mg_type(const bool & mg_type_in ) 
  //   { 
  //       _mg_type = mg_type_in; 
  //       PC pc; 
  //       KSPGetPC(this->ksp, &pc);
  //       PCMGSetType(pc, PC_MG_MULTIPLICATIVE);  // todo: change
  //   }; 

  //     SizeType  pre_smoothing_steps() const         
  //     { 
  //       return _pre_smoothing_steps; 
  //     } 

  //     SizeType  post_smoothing_steps() const        
  //     { 
  //       return _post_smoothing_steps; 
  //     }

  //     virtual bool solve(const Vector &rhs, Vector & x_0)
  //     {
  //       KSPSetUp(this->ksp);
  //       KSPSolve(this->ksp, raw_type(rhs), raw_type(x_0));
  //       KSPDestroy(&this->ksp);
  //       return true; 
  //     }

  //     virtual bool apply(const Vector &rhs, Vector & x_0) override
  //     {
  //         solve(rhs, x_0);
  //         KSPDestroy(&this->ksp);
  //         return true; 
  //     }

  //     virtual void update(const std::shared_ptr<const Matrix> &op) override
  //     {
  //         KSPSolver<Matrix, Vector, PETSC>::update(op);
  //         // TODO:: this is not very efficient, isnce we init P/R  every thime, we do LS in NLS ...
  //         this->galerkin_assembly(*op); 
  //     }

  //     virtual bool galerkin_assembly(const Matrix & A) override
  //     {
  //         PC pc;
  //         MPI_Comm            comm; 
          
  //         PetscObjectGetComm((PetscObject)raw_type(A), &comm);
  //         KSPCreate(comm, &this->ksp);
  //         KSPSetOperators(this->ksp, raw_type(A), raw_type(A));

  //         KSPSetType(this->ksp, KSPRICHARDSON); 
  //         // KSPSetType(this->ksp, KSPBCGS); 
  //         KSPSetFromOptions(this->ksp);

  //         KSPGetPC(this->ksp, &pc);
  //         PCSetType(pc, PCMG); 
  //         PCMGSetLevels(pc, this->_transfers.size() + 1, NULL);
  //         PCMGSetGalerkin(pc, PETSC_TRUE);


  //         for (SizeType k=1; k<this->_transfers.size() + 1; k++) 
  //         {
  //           Matrix P = this->_transfers[k-1].I(); 
  //           PCMGSetInterpolation(pc, k, raw_type(P));
  //         }

  //       return true; 
  //     }

  // private:
  //   std::string preconditioner_;
  //   SizeType                            _pre_smoothing_steps; 
  //   SizeType                            _post_smoothing_steps; 
  //   SizeType                            _mg_type; 
  // };

}
//////////////////////////////////////////////////////////////////////////////////////////////////

#endif //UTOPIA_PETSC_KSP_SOLVERS_HPP
