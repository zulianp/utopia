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


        // static_assert(Traits<Matrix>::Backend == utopia::PETSC_EXPERIMENTAL, "only works with petsc types");


        SNESSolver( const std::shared_ptr <LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(),
                    const Parameters params                       = Parameters() ):
                    NonLinearSolver<Matrix, Vector>(linear_solver, params)
        {

        }


        virtual ~SNESSolver()
        { 

        }




        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
           using namespace utopia;

           std::cout<<"I am inside of solve function .......... \n"; 


           SNES           snes;
           KSP            ksp; 
           PC             pc; 


           // residual 
           Vector residual = local_zeros(local_size(x)); 


            SNESCreate(PETSC_COMM_WORLD,&snes);


            // enrgy 
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

                                  // disp(x_ut); 

                                  // this is horrible copying of mats - should be fixed 
                                  Matrix jac_ut; // = sparse_mref(jac);  // because solver test is using dense matrix ... 

                                  fun->hessian(x_ut, jac_ut); 


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


            // // needs to be changed 
            // // also take into account utopia LS and PCs 
            SNESGetKSP(snes,&ksp);
            KSPSetType(ksp, KSPGMRES); 

            KSPGetPC(ksp,&pc);
            PCSetType(pc,PCLU);

            KSPSetTolerances(ksp,1e-11, 1e-11, PETSC_DEFAULT, 50000); 
            SNESSetTolerances(snes, 1e-11, 1e-11, 1e-11, 1000, 500000); 
            
            // lets also print other stuff, like energy and so on...
             SNESMonitorSet(
                 snes,
                 [](SNES snes, PetscInt iter, PetscReal res, void*) -> PetscErrorCode {
                     PrintInfo::print_iter_status({static_cast<PetscReal>(iter), res}); 
                     return 0;
                 },
                 nullptr,
                 nullptr);


            SNESSetFromOptions(snes);
            KSPSetFromOptions(ksp);


            SNESSetType(snes, SNESNEWTONTR);
         
            SNESSolve(snes, NULL, raw_type(x));
            
            PetscInt nonl_its; 
            SNESGetIterationNumber(snes, &nonl_its);


            SNESDestroy(&snes);

            std::cout<<"number of nonlinear iterations: "<< nonl_its << " \n"; 

           return true; 
       }





//TO BE DONE:
// - ksp utopia
// - params
// - set type 
// - preconditioner 
// - convert functions 
// - nicer represnetation of printout
// - exit solver
// - allocation of hessian 
//  - proper types inside of fun... 
// - utopia_DM ??? yes/no?? maybe?? 
// - nullspaces ?? near-nullspaces?? - maybe just libmesh based people ...
// 




    };
    
}

#endif // UTOPIA_SNES_HPP  
