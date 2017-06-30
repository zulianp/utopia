/*
* @Author: Alena Kopanicakova
* @Date:   2017-05-22
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-06-30
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

        // We separate cases with 1 and 2 constraints in order to avoid usless computations in single constraint case
        bool solve(const Matrix &A, const Vector &b, Vector &x)  override
        {
            if( constraints_.has_upper_bound() && constraints_.has_lower_bound())
                box_solve(A, b, x); 
            else if((constraints_.has_upper_bound() && !constraints_.has_lower_bound()) || (!constraints_.has_upper_bound() && constraints_.has_lower_bound()))
                single_bound_solve(A, b, x); 
            else
                std::cout<<"if u do not have constraints, use something else..... \n"; 
            return true; 
        }


        virtual void set_parameters(const Parameters params) override
        {
            IterativeSolver<Matrix, Vector>::set_parameters(params);
            c_ = 1.0; 
        }


        virtual bool set_box_constraints(const BoxConstraints & box)
        {
          constraints_ = box; 
          return true; 
        }

        virtual bool  get_box_constraints(BoxConstraints & box)
        {
            box = constraints_; 
            return true; 
        }


        virtual void set_c(const Scalar & c)
        {
            c_ = c; 
        }

        virtual void get_c(Scalar & c) const 
        {
            c = c_; 
        }


private: 

        bool single_bound_solve(const Matrix &A, const Vector &b, Vector &x)
        {
            using namespace utopia;

            SizeType n = local_size(A).get(0), m = local_size(A).get(1); 

            Scalar x_diff_norm;


            SizeType it = 0;
            bool converged = false; 

            Vector bound; 

            if(constraints_.has_upper_bound())
                bound = *constraints_.upper_bound(); 
            else
                bound = *constraints_.lower_bound(); 

            Vector lambda    = local_zeros(n);
            Vector active    = local_zeros(n);
            Vector G_inv     = local_zeros(n);

            Vector x_old = x;

            Vector d, prev_active;

            // active/inactive constraints 
            Matrix A_s, Ic;

            // G can be changed to something else
            Matrix G = identity(n, n);

            Matrix H;

            linear_solver_->solve(G, bound, G_inv); 
    

            if(this->verbose())
                this->init_solver("SEMISMOOTH NEWTON METHOD", {" it. ", "      || x_k - x_{k-1} ||"});
            
            while(!converged) 
            {
                d = lambda + c_ * (G * x - bound);

                if(is_sparse<Matrix>::value) 
                {
                    Ic   = local_sparse(n, m, 1); 
                    Vector o = local_values(n, 1.0); 
                    Ic = diag(o); 
                } 
                else 
                {
                    Ic   = local_values(n, m, 1.0); 
                }
               
                active   = local_zeros(n);

                {
                    Write<Vector> wa(active);
                    Read<Vector> rdp(d);                                   

                    if(constraints_.has_upper_bound())
                    {
                        Range rr = range(d);
                        for (SizeType i = rr.begin(); i != rr.end(); i++) 
                        {
                            if (d.get(i) >= 0) 
                                active.set(i, 1.0);
                        }
                    }
                    else
                    {
                        Range rr = range(d);
                        for (SizeType i = rr.begin(); i != rr.end(); i++) 
                        {
                            if (d.get(i) <= 0) 
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
                        if(this->verbose())
                            std::cout<<"SemismoothNewton:: set is not changing ... \n"; 
                        converged = true; 
                    }
                }

                prev_active = active;        

                H = A_s + Ic * A;         

                Vector rhs =  Ic * b + A_s * G_inv; 

                linear_solver_->solve(H, rhs, x); 
    
                lambda = A_s * (b - A * x);

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

        // TODO: make this solve more efficient
        bool box_solve(const Matrix &A, const Vector &b, Vector &x)  
        {
            using namespace utopia;

            SizeType n = local_size(A).get(0), m = local_size(A).get(1); 
            Scalar x_diff_norm;

            SizeType it = 0;
            bool converged = false; 

            Vector lb = *constraints_.lower_bound(); 
            Vector ub = *constraints_.upper_bound(); 

            Vector lambda_p     = local_zeros(n);
            Vector lambda_m     = local_zeros(n);
            Vector active       = local_zeros(n);
            Vector G_inv_ub     = local_zeros(n);
            Vector G_inv_lb     = local_zeros(n);
            

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


                if(is_sparse<Matrix>::value) 
                {
                    Ac_p = local_sparse(n, m, 1); 
                    Ac_m = local_sparse(n, m, 1); 

                    Ac_p  = diag(Vector(local_zeros(n))); 
                    Ac_m  = diag(Vector(local_zeros(n))); 

                    Ic   = local_identity(n, m); 
                    Vector o = local_values(n, 1.0); 
                    Ic = diag(o); 
                } 
                else 
                {
                    Ac_p = local_zeros(n); 
                    Ac_m = local_zeros(n);
                    Ic   = local_identity(n, m); 
                }


                while(!converged) 
                {
                    d_p = lambda_p + c_ * (G * x - ub);
                    d_m = lambda_m + c_ * (G * x - lb);

                    {
                        Write<Matrix> wAC_p(Ac_p);
                        Write<Matrix> wAC_m(Ac_m);
                        Write<Vector> wa(active);

                        Read<Vector> rdp(d_p);                    
                        Read<Vector> rdm(d_m);               
                        Read<Vector> rb(b);               
                        

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
                            else if (lambda_p.get(i) < 0) 
                            {
                                Ac_p.set(i, i, 0.0);
                                active.set(i, 0.0);
                            }
                            else if(lambda_m.get(i) < 0)
                            {
                                Ac_m.set(i, i, 0.0);
                                active.set(i, 0.0);
                            }

                        }
                    }


                    A_s  = diag(active); 
                    Ic   = local_identity(n, m); 
                    Ic -= A_s; 


                    if (it > 1) 
                    {
                        // active set doesn't change anymore => converged 
                        SizeType sum = norm1(prev_active - active);
                        if (sum == 0) 
                        {
                            if(this->verbose())
                                std::cout<<"SemismoothNewton:: set is not changing ... \n"; 
                            converged = true; 
                        }
                
                    }

                    prev_active = active;        
                    H = A_s + Ic * A;         

                    Vector rhs =  Ic * b + Ac_p * G_inv_ub + Ac_m * G_inv_lb; 
                    linear_solver_->solve(H, rhs, x); 
        
                    lambda_p = Ac_p * (b - A * x);
                    lambda_m = Ac_m * (A * x - b);


                    // check for change in iterates 
                    x_diff_norm = (it > 1) ? norm2(x - x_old): 9e9; 

                    
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



        std::shared_ptr <Solver>        linear_solver_;  
        BoxConstraints                  constraints_; 
        Scalar c_; 
    };

}
#endif //UTOPIA_SOLVER_SEMISMOOTH_NEWTON_HPP
