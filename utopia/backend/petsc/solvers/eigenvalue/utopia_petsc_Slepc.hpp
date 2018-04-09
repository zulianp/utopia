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


        EigenvelueProblemSlover(    const Parameters params = Parameters(), 
                                    const std::vector<std::string> problem_types    = {"hermitian", "non_hermitian", "generalized_hermitian", "generalized_non_hermitian", "generalized_hermitian_SPD_B", "generalized_hermitian_indefinite"}, 
                                    const std::vector<std::string> portions_of_spectrum    = {"largest_magnitude", "smallest_magnitude", "largest_real", "smallest_real", "largest_imaginary", "smallest_imaginary", "closest_to_target", "closest_to_target_real", "closest_to_target_imaginary", "all_in_region"},
                                    const std::vector<std::string> solver_types    = {"krylovschur", "power", "subspace", "arnoldi", "lanczos", "gd", "jd", "rqcg", "lobpcg", "ciss", "lapack", "arpack", "blzpack", "trlan", "blopex", "primme", "feast"}): 
                                    initialized_(false), 
                                    solved_(false), 
                                    number_of_eigenvalues_(1),
                                    problem_types_(problem_types), 
                                    portions_of_spectrum_(portions_of_spectrum), 
                                    solver_types_(solver_types), 
                                    max_it_(1000), 
                                    tol_(1e-12)

        {
            problem_type_           = problem_types_.at(0); 
            portion_of_spectrum_    = portions_of_spectrum_.at(0);
            solver_type_            = solver_types_.at(0); 


            verbose_                = params.verbose();
        }   


        virtual ~EigenvelueProblemSlover()
        { 
            if (initialized_)
                EPSDestroy(&eps_);
        }


        virtual void problem_type(const std::string & type)
        {
          problem_type_ = in_array(type, problem_types_) ? type : problem_types_.at(0);
        }


        virtual const std::string & problem_type() const 
        {
          return problem_type_; 
        }


        virtual void solver_type(const std::string & type)
        {
          solver_type_ = in_array(type, solver_types_) ? type : solver_types_.at(0);
        }


        virtual const std::string & solver_type() const 
        {
          return solver_type_; 
        }        


        virtual void portion_of_spectrum(const std::string & type)
        {
          portion_of_spectrum_ = in_array(type, portions_of_spectrum_) ? type : portions_of_spectrum_.at(0);
        }


        virtual const std::string &  portion_of_spectrum() const 
        {
          return portion_of_spectrum_; 
        }


        virtual void number_of_eigenvalues(const SizeType & number_of_eigenvalues)
        {
            number_of_eigenvalues_ = number_of_eigenvalues; 
        }

        virtual const SizeType & number_of_eigenvalues() const 
        {
            return number_of_eigenvalues_; 
        }

        virtual void max_it(const SizeType & max_it)
        {
            max_it_ = max_it; 
        }

        virtual const SizeType & max_it() const 
        {
            return max_it_; 
        }

        virtual void tol(const Scalar & tol)
        {
            tol_ = tol; 
        }

        virtual const Scalar & tol() const 
        {
            return tol_; 
        }


        virtual const bool & verbose() const 
        {
            return verbose_; 
        }


        virtual void verbose(const bool & verbose)  
        {
            verbose_ = verbose; 
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

            // PetscInt its; 
            // EPSGetIterationNumber(eps_, &its);
            // std::cout<<"it:  "<< its << "  \n"; 



            if(convergence_reason > 0)
                solved_ = true; 
            else
                std::cout<<"EigenSolver did not converge... \n"; 

            return true; 
        }


        virtual bool print_eigenpairs()
        {

            if(!verbose())
                return false; 

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

            if(problem_type_ == "largest_magnitude")
                EPSSetProblemType(eps_, EPS_HEP);
            else if(problem_type_ == "generalized_hermitian")
                EPSSetProblemType(eps_, EPS_GHEP);
            else if(problem_type_ == "non_hermitian")
                EPSSetProblemType(eps_, EPS_NHEP);                            
            else if(problem_type_ == "generalized_non_hermitian")
                EPSSetProblemType(eps_, EPS_GNHEP);
            else if(problem_type_ == "generalized_hermitian_SPD_B")
                EPSSetProblemType(eps_, EPS_PGNHEP);
            else
                EPSSetProblemType(eps_, EPS_GHIEP);            


            if(portion_of_spectrum_ == "largest_magnitude")
                EPSSetWhichEigenpairs(eps_, EPS_LARGEST_MAGNITUDE );
            else if(portion_of_spectrum_ == "smallest_magnitude")
                EPSSetWhichEigenpairs(eps_, EPS_SMALLEST_MAGNITUDE   );
            else if(portion_of_spectrum_ == "largest_real")
                EPSSetWhichEigenpairs(eps_,  EPS_LARGEST_REAL );
            else if(portion_of_spectrum_ == "smallest_real")
                EPSSetWhichEigenpairs(eps_, EPS_SMALLEST_REAL );
            else if(portion_of_spectrum_ == "largest_imaginary")
                EPSSetWhichEigenpairs(eps_, EPS_LARGEST_IMAGINARY );
            else if(portion_of_spectrum_ == "smallest_imaginary")
                EPSSetWhichEigenpairs(eps_, EPS_SMALLEST_IMAGINARY );
            else if(portion_of_spectrum_ == "closest_to_target")
                EPSSetWhichEigenpairs(eps_, EPS_TARGET_MAGNITUDE );            
            else if(portion_of_spectrum_ == "closest_to_target_real")
                EPSSetWhichEigenpairs(eps_, EPS_TARGET_REAL );            
            else if(portion_of_spectrum_ == "closest_to_target_imaginary")
                EPSSetWhichEigenpairs(eps_, EPS_TARGET_IMAGINARY );            
            else if(portion_of_spectrum_ == "all_in_region")
                EPSSetWhichEigenpairs(eps_, EPS_ALL );                                                        


            EPSSetType(eps_, solver_type_.c_str()); 
            EPSSetTolerances(eps_, tol_, max_it_); 

            // could be done more sophisticated... 
            EPSSetDimensions(eps_, number_of_eigenvalues_, PETSC_DEFAULT,  PETSC_DEFAULT); 


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

        SizeType number_of_eigenvalues_; 
        const std::vector<std::string> problem_types_; 
        const std::vector<std::string> portions_of_spectrum_; 
        const std::vector<std::string> solver_types_; 


        std::string portion_of_spectrum_;
        std::string eps_type_; 
        std::string problem_type_; 
        std::string solver_type_; 


        SizeType max_it_; 
        Scalar tol_; 

        bool verbose_; 


    };
    
}



#endif //UTOPIA_PETSC_SLEPC_H