/*
* @Author: Alena Kopanicakova
* @Date:   2018-03-3
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2018-03-3
*/
#ifndef UTOPIA_SNES_HPP
#define UTOPIA_SNES_HPP  

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_petsc.hpp"


#include "utopia_petsc_UTOPIA_KSP_Solver.hpp"


#include <algorithm>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petscsnes.h>

namespace utopia 
{


    // this should be potentially nonlinear solver .... 
    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
    class SNESSolver {};



    template<typename Matrix, typename Vector>
    class SNESSolver<Matrix, Vector, PETSC_EXPERIMENTAL> : virtual public NonLinearSolver<Matrix, Vector>
    {

    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef typename NonLinearSolver<Matrix, Vector>::Solver LinearSolver;

        typedef NonLinearSolver<Matrix, Vector> NonLinearSolver;


        SNESSolver( const std::shared_ptr <LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(),
                    const Parameters params = Parameters(), 
                    const std::vector<std::string> snes_types    = {"newtonls", "newtontr", "nrichardson", "ksponly", "vinewtonrsls", "vinewtonssls", "ngmres", "qn", "shell", "ngs", "ncg", "fas", "ms", "anderson"}):
                    NonLinearSolver(linear_solver, params), 
                    SNES_types(snes_types)
        {
          SNES_type_       = SNES_types.at(0); 
          set_parameters(params); 
        }


        virtual ~SNESSolver()
        { 

        }


    // TODO:: 
    virtual void set_parameters(const Parameters params) override
    {
      NonLinearSolver::set_parameters(params); 
    }


    virtual void set_snes_type(const std::string & type)
    {
      SNES_type_ = in_array(type, SNES_types) ? type : SNES_types.at(0);
    }


    virtual void set_snes_options(SNES & snes)
    {
        PetscErrorCode ierr;

        
        SNESSetFromOptions(snes);

      
        if(this->verbose())
           SNESMonitorSet(
               snes,
               [](SNES snes, PetscInt iter, PetscReal res, void*) -> PetscErrorCode {
                   PrintInfo::print_iter_status({static_cast<PetscReal>(iter), res}); 
                   return 0;
               },
               nullptr,
               nullptr);
        

        ierr = SNESSetType(snes, SNES_type_.c_str());
        ierr = SNESSetTolerances(snes, NonLinearSolver::atol(), NonLinearSolver::rtol(), NonLinearSolver::stol(), NonLinearSolver::max_it(), PETSC_DEFAULT); 
    }



    virtual void set_ksp(SNES & snes)
    {
        KSP            ksp; 
        SNESGetKSP(snes,&ksp);
        if (dynamic_cast<KSPSolver<Matrix, Vector>*>(this->linear_solver_.get()) != nullptr)
        {
          auto utopia_ksp = dynamic_cast<KSPSolver<Matrix, Vector> *>(this->linear_solver_.get()); 
          utopia_ksp->set_ksp_options(ksp); 
          utopia_ksp->attach_preconditioner(ksp); 
        }
        else
        {
          // check if our options overwrite this 
          KSPSetFromOptions(ksp); 
          KSPSetType(ksp, KSPUTOPIA); 
              
          KSPSetSolveRoutine_UTOPIA(ksp, get_ksp_solve_routine()); 
          KSPSetTolerances_UTOPIA(ksp, get_ksp_tol_routine()); 
          KSPSetGetConvergenceReason_UTOPIA(ksp, get_ksp_convergence_routine()); 


          {
            PC pc; 
            KSPGetPC(ksp, &pc);

            PCType pc_type; 
            PCGetType(pc, &pc_type); 

            // - just for the moment ... 
            // - TO BE DONE WITH UTOPIA - preconditiner ... 
            PCSetType(pc, "none");
          }


          KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
        }
    }


      bool solve(Function<Matrix, Vector> &fun, Vector &x) override
      {
        using namespace utopia;

        SNES            snes;
        MPI_Comm        comm;

        PetscObjectGetComm((PetscObject)raw_type(x), &comm);

        
        std::string method = "SNES_INTERFACE"; 
        this->init_solver(method, {"it", "|g|"}); 


        // residual 
        Vector residual = local_zeros(local_size(x));  
      
        SNESCreate(comm, &snes);

        // energy 
        SNESSetObjective(   snes, 
                              // FormObjective, 
                              [](SNES /*snes*/, Vec x, PetscReal * energy, void * ctx) -> PetscErrorCode 
                              {
                                  Function<Matrix, Vector> * fun = static_cast<Function<Matrix, Vector> *>(ctx);
                                  Vector x_ut;

                                  utopia::convert(x, x_ut); 

                                  fun->value(x_ut, *energy);

                               return 0;
                              },

                              &fun);



        // gradient 
        SNESSetFunction(    snes, 
                              raw_type(residual), 
                              // FormGradient, 
                              [](SNES snes, Vec x, Vec res, void *ctx)-> PetscErrorCode 
                              {
                                  Function<Matrix, Vector> * fun = static_cast<Function<Matrix, Vector> *>(ctx);
                                    
                                  Vector x_ut, res_ut; 
                                  utopia::convert(x, x_ut); 

                                  fun->gradient(x_ut, res_ut); 
                                  utopia::convert(res_ut, res);
                            
                                  return 0;
                              },
                              &fun);



        // hessian 
        SNESSetJacobian(    snes, 
                              snes->jacobian, 
                              snes->jacobian_pre, 
                              // FormHessian, 
                              [](SNES snes, Vec x, Mat jac, Mat prec, void *ctx)-> PetscErrorCode 
                              {
                                Function<Matrix, Vector> * fun = static_cast<Function<Matrix, Vector> *>(ctx);

                                Vector x_ut;
                                utopia::convert(x, x_ut); 

                                // this is horrible copying of mats - should be fixed 
                                Matrix jac_ut; // = sparse_mref(snes->jacobian);  


                                fun->hessian(x_ut, jac_ut); 


                                // STUPID 
                                MatDuplicate(raw_type(jac_ut), MAT_COPY_VALUES,  &jac); 
                                MatCopy(raw_type(jac_ut), jac, SAME_NONZERO_PATTERN ); 

                                MatDuplicate(jac, MAT_COPY_VALUES,  &prec); 
                                MatCopy(jac, prec, SAME_NONZERO_PATTERN ); 


                                  //  interesting....
                                  snes->jacobian = jac; 
                                  snes->jacobian_pre = prec; 

                                return 0;
                              },
                              &fun);


          set_snes_options(snes); 
          set_ksp(snes); 
       
          SNESSolve(snes, NULL, raw_type(x));

          // exit solver 
          PetscInt nonl_its; 
          SNESGetIterationNumber(snes, &nonl_its);

          SNESConvergedReason reason; 
          SNESGetConvergedReason(snes, &reason); 

          this->exit_solver(nonl_its, reason); 

          // cleaning
          SNESDestroy(&snes);

         return true; 
     }


     protected: 
      std::string SNES_type_;                                  /*!< Choice of snes types. */  
      const std::vector<std::string> SNES_types;              /*!< Valid options for SNES solver types. */  




    private: 

      std::function< void(PetscInt &,  KSPConvergedReason&) > get_ksp_convergence_routine()
      {
        std::function< void(PetscInt &,  KSPConvergedReason&) > fun = [this](PetscInt & max_it,  KSPConvergedReason & reason)
        {
          if (dynamic_cast<IterativeSolver<Matrix, Vector>*>(this->linear_solver_.get()) != nullptr)
          {
            auto ls = dynamic_cast<IterativeSolver<Matrix, Vector> *>(this->linear_solver_.get()); 
            max_it = ls->get_num_it(); 

            switch (ls->get_convergence_reason() )
            {
                // sucess
                case ConvergenceReason::CONVERGED_FNORM_ABS:
                    reason = KSP_CONVERGED_ATOL;
                    break;
                    
                case ConvergenceReason::CONVERGED_FNORM_RELATIVE:
                    reason = KSP_CONVERGED_RTOL;
                    break;
                    
                case ConvergenceReason::CONVERGED_SNORM_RELATIVE:
                    reason = KSP_CONVERGED_STEP_LENGTH;
                    break;
                        
                // fail
                case ConvergenceReason::DIVERGED_MAX_IT :
                    reason = KSP_DIVERGED_ITS;
                    break;
                    
                default :
                    reason = KSP_CONVERGED_RTOL_NORMAL;
            }
          }
          else
          {
            std::cout<<"get convergence reason is not configured for direct solvers yet... \n"; 
          }
        }; 
        return fun; 
      }


    std::function< void(const PetscReal &, const PetscReal &,  const PetscReal &,  const PetscInt &) > get_ksp_tol_routine()
    {
      std::function< void(const PetscReal &, const PetscReal &, 
                          const PetscReal &, const PetscInt &) > fun = [this](const PetscReal & rtol, 
                                                                              const PetscReal & abstol, 
                                                                              const PetscReal & dtol, 
                                                                              const PetscInt & maxits)
      {
        if (dynamic_cast<IterativeSolver<Matrix, Vector>*>(this->linear_solver_.get()) != nullptr)
        {
          auto ls = dynamic_cast<IterativeSolver<Matrix, Vector> *>(this->linear_solver_.get()); 
          ls->atol(abstol); 
          ls->rtol(rtol); 
          ls->max_it(maxits); 
        }
        else
        {
          std::cout<<"set tol is not configured for direct solvers yet... \n"; 
        }
      };
      return fun; 
    }


    std::function<void(const Mat &, const Mat &, const Vec &, Vec &)> get_ksp_solve_routine()
    {
      std::function<void(const Mat &, const Mat &, const Vec &, Vec &)> fun = [this](const Mat &A, const Mat & P, const Vec &b, Vec & x)
      {

        Vector  b_ut, x_ut; 

        // we need to get some better way how to do this
        convert(x, x_ut); 
        convert(b, b_ut); 

        const Matrix A_ut = sparse_mref(A); 

        // PreconditionedSolver - maybe in future ... 
        // this->linear_solver_->solve(A_ut, P_ut, b_ut, x_ut); 

        this->linear_solver_->solve(A_ut, b_ut, x_ut); 

        convert(x_ut, x); 
      };
      return fun; 
    }




  //TO BE DONE:
  // - prepare jacobian global mat - preallocation ... 
  // - convert functions could be more efficient ...
  // - allocation of hessian 


    };
}

#endif // UTOPIA_SNES_HPP  
