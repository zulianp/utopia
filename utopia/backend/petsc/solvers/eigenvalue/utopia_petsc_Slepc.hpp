/*
* @Author: Alena Kopanicakova
* @Date:   2018-04-06
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2018-04-06
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


    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
    class EigenvelueProblemSlover; 

    template<typename Matrix, typename Vector>
    class EigenvelueProblemSlover<Matrix, Vector, PETSC_EXPERIMENTAL>
    {

    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


        EigenvelueProblemSlover( ): 
                                    initialized_(false), 
                                    solved_(false)

        {

        }


        virtual ~EigenvelueProblemSlover()
        { 
            if (initialized_)
                EPSDestroy(&eps_);
        }


        virtual bool solve(const Matrix & A)
        {
            MPI_Comm            comm; 
            PetscObjectGetComm((PetscObject)raw_type(A), &comm);

            if (initialized_)
                reinitialize(comm);
            else
                initialize(comm); 

            EPSSetOperators(eps_, raw_type(A), NULL);
            EPSSolve(eps_); 

            EPSConvergedReason convergence_reason; 
            EPSGetConvergedReason(eps_, &convergence_reason); 

            if(convergence_reason > 0)
                solved_ = true; 
            else
                std::cout<<"Convergence was not successful... \n"; 

            return true; 
        }


        virtual bool print_eigenpairs()
        {
            if(!solved_)
            {
                std::cerr<<"You need to solve before query eigenpairs.... \n"; 
                return false; 
            }

            PetscInt nconv; 

            PetscScalar         kr,ki;
            Vector              xr,xi;

            Mat A; 

            EPSGetOperators(eps_, &A, NULL); 

            MatCreateVecs(A, NULL, & raw_type(xr));
            MatCreateVecs(A, NULL, & raw_type(xi));


            EPSGetConverged(eps_, &nconv);
            std::string method = "EIGEN_PAIRS  \n Number of converged eigenpairs: " + std::to_string(nconv); 
            PrintInfo::print_init(method, {"index", "real eigenvalue", "im eigenvalue"}); 

            if (nconv>0) 
            {
              for (auto i=0; i<nconv; i++) 
              {
                EPSGetEigenpair(eps_, i, &kr, &ki ,raw_type(xr), raw_type(xi));
                PrintInfo::print_iter_status(i, {kr, ki});
              }
              PetscPrintf(PETSC_COMM_WORLD,"\n");
            }

            return true; 
        }


        virtual void get_eigenpairs(const SizeType & i, Scalar & iegr, Scalar & eigi, Vector & vr, Vector & vi)
        {
            if(!solved_)
            {
                std::cerr<<"You need to solve before query eigenpairs.... \n"; 
                return; 
            }

            PetscInt nconv; 

            Mat A; 

            EPSGetOperators(eps_, &A, NULL); 

            if(!empty(vr))
                VecDestroy(&raw_type(vr)); 
            if(!empty(vi))
                VecDestroy(&raw_type(vi));

            MatCreateVecs(A, NULL, & raw_type(vr));
            MatCreateVecs(A, NULL, & raw_type(vi));

            vr.implementation().set_initialized(true); 
            vi.implementation().set_initialized(true); 

            EPSGetConverged(eps_, &nconv);


            if (i < nconv) 
                EPSGetEigenpair(eps_, i, &iegr, &eigi ,raw_type(vr), raw_type(vi));

        }



        virtual void get_real_eigenpair(const SizeType & i, Scalar & iegr, Vector & vr)
        {
            if(!solved_)
            {
                std::cerr<<"You need to solve before query eigenpairs.... \n"; 
                return; 
            }

            PetscInt nconv; 

            Mat A; 
            EPSGetOperators(eps_, &A, NULL); 

            if(!empty(vr))
                VecDestroy(&raw_type(vr)); 


            MatCreateVecs(A, NULL, & raw_type(vr));
            vr.implementation().set_initialized(true); 

            EPSGetConverged(eps_, &nconv);


            if (i < nconv) 
                EPSGetEigenpair(eps_, i, &iegr, NULL ,raw_type(vr), NULL);

        }



    private: 
        void initialize(const MPI_Comm & comm)
        {

            EPSCreate(comm, &eps_);
            EPSSetProblemType(eps_, EPS_HEP);

            EPSSetFromOptions(eps_);

            initialized_ = true; 
        }

        void reinitialize(const MPI_Comm & comm)
        {
            EPSDestroy(&eps_);
            initialize(comm); 

            solved_ = false; 
        }


    private: 
        EPS eps_; 
        bool initialized_; 
        bool solved_; 

        std::string eps_type_; 
        std::string problem_type_; 


    };
    
}



#endif //UTOPIA_PETSC_SLEPC_H