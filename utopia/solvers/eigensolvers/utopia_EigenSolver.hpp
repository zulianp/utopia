#ifndef UTOPIA_EIGEN_SOLVER_HPP
#define UTOPIA_EIGEN_SOLVER_HPP

#include "utopia_Preconditioner.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"

namespace utopia 
{

    template<typename Matrix, typename Vector>
    class EigenSolver: public virtual Clonable 
    {

    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


        EigenSolver(    const Parameters params = Parameters() ): 
                        number_of_eigenvalues_(1),
                        max_it_(1000), 
                        tol_(1e-12)

        {
            verbose_                = params.verbose();
        }   


        virtual ~EigenSolver()
        { 

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


        virtual bool solve(const Matrix & A) = 0; 
        virtual bool solve(const Matrix & A, const Matrix & B) = 0; 

        virtual bool print_eigenpairs() = 0; 
        virtual void get_eigenpairs(const SizeType & i, Scalar & iegr, Scalar & eigi, Vector & vr, Vector & vi) = 0; 
        virtual void get_real_eigenpair(const SizeType & i, Scalar & iegr, Vector & vr) = 0; 

        // virtual EigenSolver<Matrix, Vector> * clone() =0; 


    private: 
        SizeType number_of_eigenvalues_; 

        SizeType max_it_; 
        Scalar tol_; 

        bool verbose_; 


    };
    
}



#endif //UTOPIA_EIGEN_SOLVER_HPP