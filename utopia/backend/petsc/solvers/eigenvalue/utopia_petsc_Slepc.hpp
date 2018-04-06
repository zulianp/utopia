/*
* @Author: Alena Kopanicakova
* @Date:   2016-09-01
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-06-09
*/
#ifndef UTOPIA_PETSC_SLEPC_H
#define UTOPIA_PETSC_SLEPC_H

#include "utopia_Preconditioner.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_Smoother.hpp"

#include "utopia_Core.hpp"

#include <slepceps.h>

namespace utopia 
{

    template<typename Matrix, typename Vector>
    class Slepc
    {

    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


        Slepc( )

        {
                std::cout<<"-------- Slepc() --------- \n"; 
        }


        virtual ~Slepc()
        { 
        }


        virtual bool bla(const Matrix & A)
        {
            std::cout<<"------------- Slepc::bla bla bla -------- \n"; 

            EPS eps;


            PetscInt its, maxit, nev, nconv; 
            PetscReal tol; 

            PetscScalar    kr,ki;
            Vector            xr,xi;


            EPSType type   = EPSLANCZOS; 

            MPI_Comm            comm; 
            PetscObjectGetComm((PetscObject)raw_type(A), &comm);


            EPSCreate(comm, &eps);
            EPSSetOperators(eps, raw_type(A), NULL);
            EPSSetProblemType(eps, EPS_HEP);

            EPSSetFromOptions(eps);

            EPSSolve(eps); 

            EPSGetIterationNumber(eps,&its);
            PetscPrintf(comm," Number of iterations of the method: %D\n",its);


            EPSGetType(eps,&type);
            PetscPrintf(comm," Solution method: %s\n\n",type);

            EPSGetDimensions(eps,&nev,NULL,NULL);
            PetscPrintf(comm," Number of requested eigenvalues: %D\n",nev); 


            EPSGetTolerances(eps,&tol,&maxit);


            PetscPrintf(comm," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);


            MatCreateVecs(raw_type(A), NULL, & raw_type(xr));
            MatCreateVecs(raw_type(A), NULL, & raw_type(xi));


            EPSGetConverged(eps,&nconv);
            PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);

            if (nconv>0) 
            {
              for (PetscInt i=0;i<nconv;i++) 
              {
                EPSGetEigenpair(eps,i,&kr,&ki,raw_type(xr),raw_type(xi));

                disp(xr); 


                std::cout<<"kr: "<< kr << "  \n"; 

              }
              PetscPrintf(PETSC_COMM_WORLD,"\n");
            }


            // VecDestroy(&xr);
            // VecDestroy(&xi);

            EPSDestroy(&eps);


            return true; 
        }


    };
    
}



#endif //UTOPIA_PETSC_SLEPC_H