/*! \file utopia_NonlinSemismoothNewton.hpp
    Semismooth Newton method 
    Created by Alena Kopanicakova 
*/

#ifndef UTOPIA_SOLVER_SEMISMOOTH_NEWTON_HPP
#define UTOPIA_SOLVER_SEMISMOOTH_NEWTON_HPP

#include "utopia_Wrapper.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include <vector>


namespace utopia 
{
    template<class Matrix, class Vector>
    class SemismoothNewton : public NonLinearSolver<Matrix, Vector>  
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef typename NonLinearSolver<Matrix, Vector>::Solver Solver;

    public:
        //! SemismoothNewton solver 
        /*
          ! Linear solver required. 
        */
        SemismoothNewton(   const std::shared_ptr <Solver> &linear_solver   = std::shared_ptr<Solver>(),
                            const Parameters params                         = Parameters() ) : 
                            NonLinearSolver<Matrix, Vector>(linear_solver, params)
                { }

        //! solve Nonlin - Semismooth Newton
        /*!
          \param x_new  initial guess/output
          \param A - laplacian
          \param b - rhs 
          \param g - upper bound
          \param tol - tolerance = 1e-8
          \param maxIter = 100
        */             
        bool solve(Vector &x_new, Matrix &A, Vector &b, Vector &g)  
        {
            using namespace utopia;

            SizeType local_N = local_size(x_new).get(0);

            Scalar c = 1;
            SizeType iterations = 0;
            bool converged = false; 

            Vector lambda = local_zeros(local_N);
            Vector active = local_zeros(local_N);
            Vector Ginvg  = local_zeros(local_N);
            Vector x_old = x_new;

            Vector d, prev_active;

            // active/inactive constraints 
            Matrix Ac;
            Matrix Ic;

            //G can be changed to something else
            Matrix G = local_identity(local_N, local_N);

            Matrix M;

            if(!this->linear_solve(G, g, Ginvg)) {
                return false;
            }

            Scalar f_norm = 999999999999;

            this->init_solver("SEMISMOOTH NEWTON METHOD", {" it. ", "|| g ||"});
            
            while(!converged) 
            {
                d = lambda + c * (G * x_new - g);

                if(is_sparse<Matrix>::value) {
                    Ac = local_sparse(local_N, local_N, 1);
                    Ic = local_sparse(local_N, local_N, 1);
                } else {
                    Ac = local_zeros({local_N, local_N});
                    Ic = local_zeros({local_N, local_N});
                }
               
                {
                    Write<Matrix> wAC(Ac);
                    Write<Vector> wa(active);
                    Write<Matrix> wIC(Ic);

                    Range rr = row_range(Ac);
                    for (SizeType i = rr.begin(); i != rr.end(); i++) 
                    {
                        if (d.get(i) > 0) 
                        {
                            Ac.set(i, i, 1.0);
                            active.set(i, 1.0);
                        }
                        else 
                        {
                            Ic.set(i, i, 1.0);
                        }
                    }
                }

                if (iterations > 1) 
                {
                    SizeType sum = norm1(prev_active - active);

                    // active set doesn't change anymore => converged 
                    if (sum == 0) 
                    {
                        // fix this to be done in other way 
                        converged = this->check_convergence(iterations, 1e-15, 1, 1); 
                        return true;
                    }
                }

                prev_active = active;

                M = Ac + Ic * A;
                
                if(!this->linear_solve(M, (Ic * b + Ac * Ginvg), x_new)) {
                    return false;
                }

                lambda = Ac * (b - A * x_new);

                f_norm = norm2(x_new - x_old);
                
                // print iteration status on every iteration 
                if(this->verbose_)
                    PrintInfo::print_iter_status(iterations, {f_norm}); 

                converged = this->check_convergence(iterations, f_norm, 1, 1); 
                
                x_old = x_new;
                iterations++;
            }
            
            return true;
        }


        // just to be here 
        bool solve(Function<Matrix, Vector> &, Vector & ) override 
        {
            return false;
        }

        private:
    };

}
#endif //UTOPIA_SOLVER_SEMISMOOTH_NEWTON_HPP
