/*
* @Author: alenakopanicakova
* @Date:   2016-04-04
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-10-12
*/

#ifndef UTOPIA_JACOBI_HPP
#define UTOPIA_JACOBI_HPP
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"

namespace utopia 
{
    template<class Matrix, class Vector>
    class PointJacobi : public  Smoother<Matrix, Vector>, public IterativeSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::IterativeSolver<Matrix, Vector> Solver;
        typedef utopia::Smoother<Matrix, Vector> Smoother;
        


    public:
        /**
         * @brief      Very simple Jacobi solver.
         *
         * @param[in]  omega  The relaxation parameter.
         */
        PointJacobi(const Parameters params = Parameters())
        {
            set_parameters(params); 
        }
        

        bool smooth(const Matrix &A, const Vector &rhs, Vector &x) override
        {        
            for(SizeType it = 0; it < this->sweeps(); it++)
            {
                sweep(A, rhs, x); 
            }
            return true; 
        }

        bool solve(const Matrix &A, const Vector &rhs, Vector &x) override
        {
            SizeType it = 0; 
            Scalar r_norm = 9999; 
            this->init_solver("Point Jacobi", {"it. ", "||r||" }); 

            while(it++ < this->max_it() && r_norm > this->rtol())
            {
                sweep(A, rhs, x); 
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {r_norm}); 
            }
            return 0; 
        }


        void set_parameters(const Parameters params) override
        {
            Smoother::set_parameters(params); 
            Solver::set_parameters(params); 
        }

    protected:

        /**
         * @brief      This implementation could be way nicer, 
                        but it is here, just to test if ML/MG works in parallel as whole thing
         *
         * @param[in]  A     The stifness  matrix.
         * @param[in]  rhs   The rhs.
         * @param      x     The initial guess/solution.
         *
         *
         * @return     { }
         */
        bool sweep(const Matrix &A, const Vector &rhs, Vector &x)
        {
            Vector diag_A = diag(A); 
            Vector diag_A_inv = 1 / diag_A; 
            
            // prevents system from being indefinite 
            check_indef(diag_A_inv); 

            // lower and upper part of A 
            Matrix LU = A;  
            Matrix D_inv, diag_A_mat = sparse(A.size().get(0), A.size().get(1), 1); 

            diag_A_mat = diag(diag_A); 
            LU -= diag_A_mat; 
             
            D_inv = diag(diag_A_inv);

            Vector r = rhs -  (LU * x); 
            x = D_inv * r; 
            
            return 0; 
        } 


        /**
         * @brief      Checks if there is a zero in the vector, if yes turn it into 1.
         *
         * @param      diag_A  { D_{-1}}
         *
         * @return     {  }
         */
        bool check_indef(Vector &diag_A) 
        {
            {
                ReadAndWrite<Vector> w(diag_A);
                Range rr = range(diag_A);

                for (SizeType i = rr.begin(); i != rr.end(); i++) 
                {
                    if(diag_A.get(i) == 0)
                    {
                        diag_A.set(i, 1);
                    }
                }

            }
            return 0; 
        } 
    
};

}

#endif //UTOPIA_JACOBI_HPP

