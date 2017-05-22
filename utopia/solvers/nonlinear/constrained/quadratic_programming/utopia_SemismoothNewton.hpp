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
    class SemismoothNewton : public IterativeSolver<Matrix, Vector>  
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> Solver;
        typedef utopia::BoxConstraints<Vector>      BoxConstraints;

    public:

        SemismoothNewton(const std::shared_ptr <Solver> &linear_solver   = std::shared_ptr<Solver>(),
                         const Parameters params                         = Parameters() ) : 
                            linear_solver_(linear_solver)

        { 
            set_parameters(params);
        }


        bool solve(const Matrix &A, const Vector &b, Vector &x)  
        {
            using namespace utopia;

            SizeType n = local_size(A).get(0), m = local_size(A).get(1); 

            Scalar c = 1, x_diff_norm;


            SizeType it = 0;
            bool converged = false; 

            Vector lb = *constraints_->lower_bound(); 
            Vector ub = *constraints_->upper_bound(); 


            Vector lambda_p   = local_zeros(n);
            Vector lambda_m   = local_zeros(n);
            Vector active   = local_zeros(n);
            Vector G_inv_ub = local_zeros(n);
            Vector G_inv_lb = local_zeros(n);
            

            Vector x_old = x;

            Vector d_p, d_m, prev_active;

            // active/inactive constraints 
            Matrix Ac_p, Ac_m, A_s, Ic;

            // G can be changed to something else
            Matrix G = identity(n, n);

            Matrix H;

            linear_solver_->solve(G, ub, G_inv_ub); 
            linear_solver_->solve(G, lb, G_inv_lb); 
            


            if(this->verbose())
                this->init_solver("SEMISMOOTH NEWTON METHOD", {" it. ", "      || x_k - x_{k-1} ||"});
            
            while(!converged) 
            {
                d_p = lambda_p + c * (G * x - ub);
                d_m = lambda_m + c * (G * x - lb);

                if(is_sparse<Matrix>::value) 
                {
                    Ac_p = local_sparse(n, m, 1); 
                    Ac_m = local_sparse(n, m, 1); 

                    Ic   = local_sparse(n, m, 1); 
                    Vector o = local_values(n, 1.0); 
                    Ic = diag(o); 
                } 
                else 
                {
                    Ac_p = local_zeros(n); 
                    Ac_m = local_zeros(n);
                    Ic   = local_values(n, m, 1.0); 
                }
               
                active   = local_zeros(n);

                {
                    Write<Matrix> wAC_p(Ac_p);
                    Write<Matrix> wAC_m(Ac_m);
                    Write<Vector> wa(active);

                    Read<Vector> rdp(d_p);                    
                    Read<Vector> rdm(d_m);                    

                    Range rr = row_range(Ac_p);
                    for (SizeType i = rr.begin(); i != rr.end(); i++) 
                    {
                        if (d_p.get(i) >= 0) 
                        {
                            Ac_p.set(i, i, 1.0);
                            active.set(i, 1.0);
                        }
                        else if(d_m.get(i) <= 0)
                        {
                            Ac_m.set(i, i, 1.0);
                            active.set(i, 1.0);
                        }

                    }
                }

                A_s  = diag(active); 
                Ic -= A_s; 

                if (it > 1) 
                {
                    // active set doesn't change anymore => converged 
                    SizeType sum = norm1(prev_active - active);
                    if (sum == 0) 
                    {
                        std::cout<<"set is not changing ... \n"; 
                        converged = true; 
                    }
            
                }

                prev_active = active;        

                H = A_s + Ic * A;         

                Vector rhs =  Ic * b + Ac_p * G_inv_ub + Ac_m * G_inv_lb; 

                linear_solver_->solve(H, rhs, x); 
    
                lambda_p = Ac_p * (b - A * x); 
                lambda_m = Ac_m * (A * x - b);

                x_diff_norm = norm2(x - x_old);
                
                // print iteration status on every iteration 
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {x_diff_norm}); 

                if(!converged)
                    converged = this->check_convergence(it, x_diff_norm, 1, 1); 
                
                x_old = x;
                it++;
            }
            
            return true;
        }


        virtual void set_parameters(const Parameters params) override
        {
            IterativeSolver<Matrix, Vector>::set_parameters(params);
        }


        virtual bool set_box_constraints(const std::shared_ptr<BoxConstraints> & box)
        {
          constraints_ = box; 
          return true; 
        }

        virtual std::shared_ptr<BoxConstraints> get_box_constraints()
        {
          return constraints_; 
        }


        private:
            std::shared_ptr <Solver>        linear_solver_;  
            std::shared_ptr<BoxConstraints> constraints_; 
    };

}
#endif //UTOPIA_SOLVER_SEMISMOOTH_NEWTON_HPP
