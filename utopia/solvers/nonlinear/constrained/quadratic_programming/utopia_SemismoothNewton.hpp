/*
* @Author: Alena Kopanicakova
* @Date:   2017-05-22
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-03
*/


#ifndef UTOPIA_SOLVER_SEMISMOOTH_NEWTON_HPP
#define UTOPIA_SOLVER_SEMISMOOTH_NEWTON_HPP

#include "utopia_Wrapper.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include <vector>
#include <memory>

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

        bool solve(Vector &x, const Matrix &A, const Vector &b, const Vector &g)  
        {
            std::cerr << "[Warning][Deprecated] SemismoothNewton: use the new box constraint interface. This method will be removed shortly" << std::endl;
            std::cout << "[Warning][Deprecated] SemismoothNewton: use the new box constraint interface. This method will be removed shortly" << std::endl;

            set_box_constraints(make_upper_bound_constraints(std::make_shared<Vector>(g)));
            return solve(A, b, x);
        }


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
        // We separate cases with 1 and 2 constraints in order to avoid usless computations in single constraint case
        bool single_bound_solve(const Matrix &A, const Vector &b, Vector &x_new)
        {

            using namespace utopia;
            
            bool is_upper_bound = constraints_.has_upper_bound();
           
            Vector g; 
            if(is_upper_bound)
                g = *constraints_.upper_bound(); 
            else
                g = *constraints_.lower_bound(); 

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

            if(!linear_solver_->solve(G, g, Ginvg))
                return false;

            Scalar f_norm = 9e9;

            if(this->verbose())
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

                    if(is_upper_bound) {
                        for (SizeType i = rr.begin(); i != rr.end(); i++) 
                        {
                            if (d.get(i) >= 0) 
                            {
                                Ac.set(i, i, 1.0);
                                active.set(i, 1.0);
                            }
                            else 
                            {
                                Ic.set(i, i, 1.0);
                            }
                        }
                    } else { //is_lower_bound
                        for (SizeType i = rr.begin(); i != rr.end(); i++) 
                        {
                            if (d.get(i) <= 0) 
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
            
                if(!linear_solver_->solve(M, (Ic * b + Ac * Ginvg), x_new))
                    return false;

                lambda = Ac * (b - A * x_new);

                f_norm = norm2(x_new - x_old);
                
                // print iteration status on every iteration 
                if(this->verbose())
                    PrintInfo::print_iter_status(iterations, {f_norm}); 

                converged = this->check_convergence(iterations, f_norm, 1, 1); 
                
                x_old = x_new;
                iterations++;
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

            const Vector &lb = *constraints_.lower_bound(); 
            const Vector &ub = *constraints_.upper_bound(); 

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
            Matrix G = local_identity(n, n);

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
                            //Why these next lines???? @Alena
                            // else if (lambda_p.get(i) < 0) 
                            // {
                            //     Ac_p.set(i, i, 0.0);
                            //     active.set(i, 0.0);
                            // }
                            // else if(lambda_m.get(i) < 0)
                            // {
                            //     Ac_m.set(i, i, 0.0);
                            //     active.set(i, 0.0);
                            // }

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
