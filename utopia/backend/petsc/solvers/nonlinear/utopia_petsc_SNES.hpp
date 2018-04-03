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

        typedef utopia::NonLinearSolver<Matrix, Vector> NonLinearSolver;


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
          configure_KSP_utopia(this->linear_solver_, ksp); 
        }
    }


      bool solve(Function<Matrix, Vector> &fun, Vector &x) override
      {
        using namespace utopia;

        SNES            snes;
                
        std::string method = "SNES_INTERFACE"; 
        this->init_solver(method, {"it", "|g|"}); 


        if (dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun) != nullptr)
        {
          PETSCUtopiaNonlinearFunction<Matrix, Vector> * fun_petsc = dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun);
          fun_petsc->getSNES(snes); 
        }
        else
        {
          MPI_Comm        comm;
          PetscObjectGetComm((PetscObject)raw_type(x), &comm);

           // residual 
          Vector residual = local_zeros(local_size(x));  
        
          SNESCreate(comm, &snes);

          // energy 
          SNESSetObjective( snes, 
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
          SNESSetFunction( snes, 
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
          SNESSetJacobian( snes,  snes->jacobian,  snes->jacobian_pre, 
                          // FormHessian, 
                          [](SNES snes, Vec x, Mat jac, Mat prec, void *ctx)-> PetscErrorCode 
                          {
                            Function<Matrix, Vector> * fun = static_cast<Function<Matrix, Vector> *>(ctx);

                            Vector x_ut;
                            utopia::convert(x, x_ut); 

                            PetscBool  assembled1, assembled2; 
                            MatAssembled(snes->jacobian, &assembled1); 

                            Matrix jac_ut; 

                            if(assembled1)
                              jac_ut.implementation().wrap(snes->jacobian);  

                            Matrix jac_ut_prec; 
                            MatAssembled(snes->jacobian_pre, &assembled2); 

                            if(assembled2)
                              jac_ut_prec.implementation().wrap(snes->jacobian_pre); 


                            bool flg = fun->hessian(x_ut, jac_ut, jac_ut_prec); 

                            if(!flg)
                            { 
                              fun->hessian(x_ut, jac_ut); 
                              MatCopy(raw_type(jac_ut), snes->jacobian_pre, SAME_NONZERO_PATTERN); 
                            }

                            if(!assembled1)
                            {
                              MatCopy(raw_type(jac_ut), snes->jacobian, SAME_NONZERO_PATTERN );                                 
                            }

                            if(!assembled2 || flg)
                            {
                              MatCopy(raw_type(jac_ut), snes->jacobian_pre, SAME_NONZERO_PATTERN );                                 
                            }

                            return 0;
                          },
                          &fun);
        }

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

  // TO BE DONE:
  // - convert functions could be more efficient ...
  // - store snes and destroy only when needed 

    };
}

#endif // UTOPIA_SNES_HPP  
