#ifndef UTOPIA_CONSTRAINT_QP_MPRGP
#define UTOPIA_CONSTRAINT_QP_MPRGP

#include <string>
#include "utopia_BoxConstraints.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_For.hpp"
#include "cuda_profiler_api.h"

namespace  utopia
{

    template<class Matrix, class Vector>
    class MPGRP final:  public MatrixFreeQPSolver<Vector>, public QPSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        using Solver  = utopia::LinearSolver<Matrix, Vector>;
        using ForLoop = utopia::ParallelFor<Traits<Vector>::Backend>;

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
                Scalar pAp, beta_beta, fi_fi, gp_dot; 

                SizeType it =0; 
                bool converged= false; 
                Scalar gnorm; 

                Scalar alpha_cg, alpha_f, beta_sc; 
	
                this->get_projection(x, *lb, *ub, Ax); 
                x = Ax; 


                // cudaProfilerStart();

                A.apply(x, Ax);
                g = Ax - rhs; 
                this->get_fi(x, g, *lb, *ub, fi); 
                this->get_beta(x, g, *lb, *ub, beta); 

                gp = fi + beta; 
                p = fi; 

                dots(beta, beta, beta_beta, fi, fi, fi_fi); 

                // cudaProfilerStop();
            

                cudaProfilerStart();

                for(auto i=0; i < 50; i++){

                    // A.apply(p, Ap);

                    // dots(   p, Ap, pAp, 
                    //         g, p, gp_dot);    


                    A.apply(p, Ap);

                    dots(   p, Ap, pAp, 
                            g, p, gp_dot);                            

                    alpha_cg = gp_dot/pAp;
                    // y = x - alpha_cg*p;      
                    alpha_f = get_alpha_f(x, p, *lb, *ub, alpha_f1, alpha_f2);        

                    this->get_fi(x, g, *lb, *ub, fi); 
                    beta_sc = dot(fi,Ap)/pAp;               

                    // Vector help = x;
                    this->get_projection(x, *lb, *ub, x);       
                    g = Ax - rhs;                 

                    // g = g - alpha_f*Ap;

                    // this->get_fi(x, g, *lb, *ub, fi); 
                    this->get_beta(x, g, *lb, *ub, beta);                     

                    gp = fi+beta;

                    dots(   beta, beta, beta_beta, 
                            fi, fi, fi_fi, 
                            gp, gp, gnorm); 
                }

                cudaProfilerStop();

                return false;


                // cudaProfilerStop();

                // return false; 

                // this->get_fi(x, g, *lb, *ub, fi); 
                // this->get_beta(x, g, *lb, *ub, beta); 
                
                // gp = fi + beta; 
                // p = fi; 

                // dots(beta, beta, beta_beta, fi, fi, fi_fi); 

                while(!converged)
                {
                    // dots(beta, beta, beta_beta, fi, fi, fi_fi); 

                    if(beta_beta <= (gamma*gamma * fi_fi))
                    {
                        A.apply(p, Ap);

                        dots(   p, Ap, pAp, 
                                g, p, gp_dot); 

                        alpha_cg = gp_dot/pAp;
                        y = x - alpha_cg*p;

                        alpha_f = get_alpha_f(x, p, *lb, *ub, alpha_f1, alpha_f2);         

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

                    dots(   beta, beta, beta_beta, 
                            fi, fi, fi_fi, 
                            gp, gp, gnorm); 
                    
                    gnorm = std::sqrt(gnorm); 
                    it++; 

                    if(this->verbose()){
                        PrintInfo::print_iter_status(it, {gnorm});
                    }

                    converged = this->check_convergence(it, gnorm, 1, 1);
                }

                // cudaProfilerStop();
    
                return true;
            }

          
        public:

            void get_fi(const Vector &x, const Vector &g, const Vector &lb, const Vector &ub, Vector & fi) const
            {
                if(empty(fi)){
                    fi = local_values(local_size(x), 0); 
                }

                {
                    auto d_lb = const_device_view(lb);
                    auto d_ub = const_device_view(ub);
                    auto d_x  = const_device_view(x);
                    auto d_g  = const_device_view(g);

                    auto d_fi = device_view(fi);
                    
                    parallel_each_write(fi, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        Scalar li = d_lb.get(i);
                        Scalar ui = d_ub.get(i);
                        Scalar xi = d_x.get(i);

                        if(li < xi && xi < ui){
                            return d_g.get(i);
                        }
                        else{
                            return 0.0;
                        }

                    });

                }
            }


            Scalar get_alpha_f(const Vector &x, const Vector &p, const Vector &lb, const Vector &ub, Vector & alpha_f1, Vector & alpha_f2) const
            {
                if(empty(alpha_f1)){
                    alpha_f1 = local_values(local_size(x), 1e15); 
                }

                if(empty(alpha_f2)){
                    alpha_f2 = local_values(local_size(x), 1e15); 
                }

                {
                    auto d_lb = const_device_view(lb);
                    auto d_ub = const_device_view(ub);
                    auto d_x  = const_device_view(x);
                    auto d_p  = const_device_view(p);

                    auto d_alpha_f1 = device_view(alpha_f1);
                    auto d_alpha_f2 = device_view(alpha_f2);
                    
                    parallel_each_write(alpha_f1, UTOPIA_LAMBDA(const SizeType i) -> Scalar
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

                    parallel_each_write(alpha_f2, UTOPIA_LAMBDA(const SizeType i) -> Scalar
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

                return multi_min(alpha_f1, alpha_f2); 
            }


            void get_beta(const Vector &x, const Vector &g, const Vector &lb, const Vector &ub, Vector & beta) const
            {
                if(empty(beta)){
                    beta = local_values(local_size(x), 0.0); 
                }

                {
                    auto d_lb = const_device_view(lb);
                    auto d_ub = const_device_view(ub);
                    auto d_x  = const_device_view(x);
                    auto d_g  = const_device_view(g);

                    // auto d_beta = device_view(beta);

                    parallel_each_write(beta, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        Scalar li = d_lb.get(i);
                        Scalar ui = d_ub.get(i);
                        Scalar xi = d_x.get(i);
                        Scalar gi = d_g.get(i);

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

        private:

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

                if(!empty(alpha_f1)) {
                    alpha_f1 = zero_expr;
                }            

                if(!empty(alpha_f2)) {
                    alpha_f2 = zero_expr;
                }                                            
            }


        private:
            Vector fi, beta, gp, p, y, Ap, Abeta, Ax, g, alpha_f1, alpha_f2; 

            Scalar eps_eig_est_; 
            SizeType power_method_max_it_; 

    };
}

#endif //UTOPIA_CONSTRAINT_QP_MPRGP
