#ifndef UTOPIA_CONSTRAINT_QP_MPRGP
#define UTOPIA_CONSTRAINT_QP_MPRGP

#include <string>
#include "utopia_BoxConstraints.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_DeviceView.hpp"

namespace  utopia
{

    template<class Matrix, class Vector>
    class MPGRP final:  public MatrixFreeQPSolver<Vector>, public QPSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                   Scalar;
        typedef utopia::LinearSolver<Matrix, Vector>    Solver;

        public:
            MPGRP(): eps_eig_est_(1e-1), power_method_max_it_(10)
            {

            }

            void read(Input &in) override
            {
                MatrixFreeQPSolver<Vector>::read(in);
                QPSolver<Matrix, Vector>::read(in);
            }


            void print_usage(std::ostream &os) const override
            {
                MatrixFreeQPSolver<Vector>::print_usage(os);
                QPSolver<Matrix, Vector>::print_usage(os);
            }



            MPGRP * clone() const override
            {
                return new MPGRP(*this);
            }

            // void set_preconditioner(const std::shared_ptr<Preconditioner<Vector> > &precond)
            // {
            //     precond_ = precond;
            // }

            bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) override
            {
                this->fill_empty_bounds(); 
                auto &box = this->get_box_constraints();


                init(local_size(rhs));
                return aux_solve(A, rhs, sol, box);
            }

            bool solve(const Matrix &A, const Vector &rhs, Vector &sol) override
            {
                auto A_op_ptr = utopia::op_ref(A);
                this->fill_empty_bounds(); 
                auto &box = this->get_box_constraints();

                init(local_size(rhs));
                return aux_solve(*A_op_ptr, rhs, sol, box);
            }


            void set_eig_comp_tol(const Scalar & eps_eig_est)
            {
                eps_eig_est_ = eps_eig_est; 
            }


        private:
            bool aux_solve(const Operator<Vector> &A, const Vector &rhs, Vector &x, const BoxConstraints<Vector> & constraints)
            {
                const auto &ub = constraints.upper_bound();
                const auto &lb = constraints.lower_bound();

                if(this->verbose()){
                    this->init_solver("MPGRP", {"it", "|| g ||"});
                }

                const Scalar gamma = 1.0; 
                const Scalar alpha_bar = 1.95/this->get_normA(A, local_size(rhs)); 
                Scalar pAp; 

                SizeType it =0; 
                bool converged= false; 
                Scalar gnorm; 

                Scalar alpha_cg, alpha_f, beta_sc; 
	
                this->get_projection(x, *lb, *ub, Ax); 
                x = Ax; 

                A.apply(x, Ax);
                g = Ax - rhs; 

                this->get_fi(x, g, *lb, *ub, fi); 
             
                this->get_beta(x, g, *lb, *ub, beta); 
                gp = fi + beta; 
                p = fi; 

                while(!converged)
                {
                    if(dot(beta, beta) <= (gamma*gamma * dot(fi,fi)))
                    {
                        A.apply(p, Ap);

                        // curvature condition check?? 
                        pAp = dot(p, Ap); 

                        // if(pAp <= 0.0)
                        // {   
                        //     if(this->verbose())
                        //     {
                        //         std::cout<<"MPRGP::curvature condition violated, terminating.... \n"; 
                        //     }
                        //     return true;
                        // }


                        alpha_cg = dot(g, p)/pAp;
                        y = x - alpha_cg*p;

                        alpha_f = get_alpha_f(x, p, *lb, *ub);         

                        if(alpha_cg <= alpha_f)
                        {
                            x = y; 
                            g = g - alpha_cg*Ap;
                            this->get_fi(x, g, *lb, *ub, fi); 
                            beta_sc = dot(fi,Ap)/pAp;
                            p = fi - beta_sc*p;                                               
                        }
                        else
                        {
                            x = x-alpha_f*p;
                            g = g - alpha_f*Ap;
                            this->get_fi(x, g, *lb, *ub, fi); 

                            Vector help = x - (alpha_bar * fi); 
                            this->get_projection(help, *lb, *ub, x);   
              
                            A.apply(x, Ax);
                            g = Ax - rhs; 
                            this->get_fi(x, g, *lb, *ub, p); 
                        }

                    }
                    else
                    {
                        A.apply(beta, Abeta);
                        alpha_cg = dot(g, beta)/dot(beta, Abeta);
                        x = x - alpha_cg*beta;
                        g = g - alpha_cg*Abeta;

                        this->get_fi(x, g, *lb, *ub, p); 
                    }

                    this->get_fi(x, g, *lb, *ub, fi); 
                    this->get_beta(x, g, *lb, *ub, beta); 

                    gp = fi+beta;
                    
                    gnorm = norm2(gp); 
                    it++; 

                    if(this->verbose()){
                        PrintInfo::print_iter_status(it, {gnorm});
                    }

                    converged = this->check_convergence(it, gnorm, 1, 1);
                }
    
                return true;
            }

          

        private:

            void get_fi(const Vector &x, const Vector &g, const Vector &lb, const Vector &ub, Vector & fi) const
            {
                fi = local_values(local_size(x), 0); 

                {
                    Read<Vector> r_ub(ub), r_lb(lb), r_x(x), r_g(g);
                    Write<Vector> w1(fi); 
                    
                    each_write(fi, [&ub, &lb, &x, &g](const SizeType i) -> double 
                    {
                        Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);
                        if(li < xi && xi < ui){
                            return g.get(i);
                        }
                        else{
                            return 0.0; 
                        }
                    }   );
                }
            }


            Scalar get_alpha_f(const Vector &x, const Vector &p, const Vector &lb, const Vector &ub) const
            {
                Vector alpha_f1 = local_values(local_size(x), 1e15); 
                Vector alpha_f2 = local_values(local_size(x), 1e15); 

                {
                    auto d_lb = const_device_view(lb);
                    auto d_ub = const_device_view(ub);
                    auto d_x  = const_device_view(x);
                    auto d_p  = const_device_view(p);

                    parallel_each_write(alpha_f1, [d_lb, d_x, d_p](const SizeType i) -> double 
                    {
                        Scalar li = d_lb.get(i); 
                        Scalar xi = d_x.get(i); 
                        Scalar pi = d_p.get(i);

                        if(pi > 0)
                        {
                            return (xi-li)/pi; 
                        }
                        else
                        {
                            return 1e15; 
                        }

                    });

                    parallel_each_write(alpha_f2, [d_ub, d_x, d_p](const SizeType i) -> double 
                    {
                        Scalar ui = d_ub.get(i); 
                        Scalar xi = d_x.get(i); 
                        Scalar pi = d_p.get(i);

                        if(pi < 0)
                        {
                            return (xi-ui)/pi; 
                        }
                        else
                        {
                            return 1e15; 
                        }

                    });

                }

                Scalar alpha_f2_min = min(alpha_f2); 
                Scalar alpha_f1_min = min(alpha_f1); 
                return std::min(alpha_f1_min, alpha_f2_min); 
            }


            void get_beta(const Vector &x, const Vector &g, const Vector &lb, const Vector &ub, Vector & beta) const
            {
                beta = local_values(local_size(x), 0.0); 

                {
                    Read<Vector> r_ub(ub), r_lb(lb), r_x(x), r_g(g);
                    Write<Vector> w1(beta); 
                    
                    each_write(beta, [&ub, &lb, &x, &g](const SizeType i) -> double 
                    {
                        Scalar li = lb.get(i); 
                        Scalar ui = ub.get(i); 
                        Scalar xi = x.get(i); 
                        Scalar gi = g.get(i);

                        if(std::abs(li -  xi) < 1e-14)
                        {
                            return std::min(0.0, gi); 
                        }
                        else if(std::abs(ui -  xi) < 1e-14)
                        {
                            return std::max(0.0, gi); 
                        }
                        else
                        {
                            return 0.0; 
                        }

                    });
                }
            }


            Scalar get_normA(const Operator<Vector> &A, const SizeType & n_loc)
            {
                // Super simple power method to estimate the biggest eigenvalue 
                Vector y_old; 
                Vector y = local_values(n_loc, 1.0); 
                SizeType it = 0; 
                bool converged = false; 
                Scalar gnorm, lambda = 0.0, lambda_old; 

                while(!converged)
                {
                    y_old = y; 
                    A.apply(y_old, y);
                    y  = Scalar(1.0/Scalar(norm2(y)))*y; 
                    gnorm = norm2(y - y_old);

                    lambda_old = lambda; 

                    A.apply(y, y_old);
                    lambda = dot(y, y_old);
                    
                    converged  = ((gnorm < eps_eig_est_) || (std::abs(lambda_old-lambda) < eps_eig_est_) || it > power_method_max_it_) ?  true: false; 
                    
                    it=it+1;  
                }

                if(this->verbose())
                    std::cout<<"Power method converged in "<< it << " iterations. Largest eig: "<< lambda << "  \n"; 

                return lambda; 
            }


            void init(const SizeType &ls)
            {
                auto zero_expr = local_zeros(ls);

                //resets all buffers in case the size has changed
                if(!empty(fi)) {
                    fi = zero_expr;
                }

                if(!empty(beta)) {
                    beta = zero_expr;
                }

                if(!empty(gp)) {
                    gp = zero_expr;
                }                

                if(!empty(p)) {
                    p = zero_expr;
                }

                if(!empty(y)) {
                    y = zero_expr;
                }

                if(!empty(Ap)) {
                    Ap = zero_expr;
                }

                if(!empty(Abeta)) {
                    Abeta = zero_expr;
                }

                if(!empty(Ax)) {
                    Ax = zero_expr;
                }            

                if(!empty(g)) {
                    g = zero_expr;
                }                            
            }


        private:
            Vector fi, beta, gp, p, y, Ap, Abeta, Ax, g; 
            Scalar eps_eig_est_; 
            SizeType power_method_max_it_; 

    };
}

#endif //UTOPIA_CONSTRAINT_QP_MPRGP
