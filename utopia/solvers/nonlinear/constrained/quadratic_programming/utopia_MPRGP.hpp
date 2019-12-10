#ifndef UTOPIA_CONSTRAINT_QP_MPRGP
#define UTOPIA_CONSTRAINT_QP_MPRGP

#include <string>
#include "utopia_BoxConstraints.hpp"
#include "utopia_QPSolver.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_For.hpp"
#include "utopia_Allocations.hpp"
#include "utopia_Algorithms.hpp"

//#include "cuda_profiler_api.h"

namespace  utopia
{

    /**
     * @brief This function is implementation of: Box constrianed quadratic programming with proportioning and projection, Dostal, 1997, SIAM, Optimization
     * @details 
     */
    template<class Matrix, class Vector>
    class MPGRP final:  public OperatorBasedQPSolver<Matrix, Vector>
    {
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;

        
        using Solver   = utopia::LinearSolver<Matrix, Vector>;
        using ForLoop  = utopia::ParallelFor<Traits<Vector>::Backend>;

        public:
            MPGRP(): eps_eig_est_(1e-1), power_method_max_it_(10), initialized_(false), loc_size_(0)
            {}

            void read(Input &in) override
            {
                OperatorBasedQPSolver<Matrix, Vector>::read(in);
                in.get("eig_comp_tol", eps_eig_est_);
                in.get("power_method_max_it", power_method_max_it_);
            }

            void print_usage(std::ostream &os) const override
            {
                OperatorBasedQPSolver<Matrix, Vector>::print_usage(os);
                this->print_param_usage(os, "eig_comp_tol", "double", "Tolerance of eigen solver.", "1e-1");
                this->print_param_usage(os, "power_method_max_it", "int", "Maximum number of iterations used inside of power method.", "10");
            }


            MPGRP * clone() const override
            {
                return new MPGRP(*this);
            }


            void update(const Operator<Vector> &A) override
            {
                SizeType loc_size_rhs = A.local_size().get(0);
                if(!initialized_ || !A.comm().conjunction(loc_size_ == loc_size_rhs)) {
                    init_memory(loc_size_rhs);
                }
            }            


            bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) override
            {
                this->fill_empty_bounds(local_size(sol));
                auto &box = this->get_box_constraints();

                this->update(A); 
                return aux_solve(A, rhs, sol, box);
            }


            void set_eig_comp_tol(const Scalar & eps_eig_est)
            {
                eps_eig_est_ = eps_eig_est;
            }


        private:
            bool aux_solve(const Operator<Vector> &A, const Vector &rhs, Vector &x, const BoxConstraints<Vector> & constraints)
            {
                // UTOPIA_NO_ALLOC_BEGIN("MPRGP");
                // //cudaProfilerStart();

                // disp(x, "x_before"); 
                // std::cout<<"x_before, norm: "<< norm_infty(x) << "  \n";

                const auto &ub = constraints.upper_bound();
                const auto &lb = constraints.lower_bound();

                // disp(*A, "A");
                // disp(rhs, "rhs");

                // disp(*lb, "lb");
                // disp(*ub, "ub");
                // exit(0);

                std::cout<<"max_in: "<< max(*ub) << "   \n"; 
                std::cout<<"min_in: "<< min(*lb) << "   \n"; 


                if(this->verbose()){
                    this->init_solver("MPGRP", {"it", "|| g ||"});
                }

                // const Scalar gamma = 1.0;
                const Scalar gamma = this->get_normA(A, local_size(rhs));
                // Scalar alpha_bar = 1.95/this->get_normA(A, local_size(rhs));
                Scalar alpha_bar = this->get_normA(A, local_size(rhs));

                Scalar pAp, beta_beta, fi_fi, gp_dot;

                SizeType it =0;
                bool converged= false;
                Scalar gnorm;

                Scalar alpha_cg, alpha_f, beta_sc;

                this->get_projection(x, *lb, *ub, Ax);
                x = Ax;

                help_f1 = *ub - *lb; 
                Scalar delta = min(help_f1) - 1e-5; 

                // disp(x, "x_projection"); 

                A.apply(x, Ax);

                // disp(Ax, "Ax"); 
                // disp(rhs, "rhs"); 


                g = Ax - rhs;

                this->get_fi(x, g, *lb, *ub, fi);
                this->get_beta(x, g, *lb, *ub, beta);

                // disp(g, "g0"); 
                // disp(beta, "beta0"); 
                // disp(x, "x"); 

                // exit(0);


                gp = fi + beta;
                p = fi;

                dots(beta, beta, beta_beta, fi, fi, fi_fi);

                // beta_beta = norm_infty(beta);

                while(!converged)
                {
                    beta_beta = norm_infty(beta);
                    beta_beta=beta_beta*beta_beta; 
                    if(beta_beta <= (gamma*gamma * fi_fi))
                    {
                        std::cout<<"-------- beta beta \n"; 
                        A.apply(p, Ap);

                        dots(p, Ap, pAp, g, p, gp_dot);

                        alpha_cg = gp_dot/pAp;
                        y = x - alpha_cg*p;
                        alpha_f = get_alpha_f(x, p, *lb, *ub, help_f1, help_f2);

                        if(alpha_cg <= alpha_f)
                        {
                            std::cout<<"--------alpha_cg <= alpha_f\n"; 
                            x = y;
                            // disp(x, "x_1"); 

                            g = g - alpha_cg*Ap;
                            this->get_fi(x, g, *lb, *ub, fi);
                            beta_sc = dot(fi,Ap)/pAp;
                            p = fi - beta_sc*p;
                        }
                        else
                        {
                            std::cout<<"--------loop 2 \n"; 
                            x = x-alpha_f*p;
                            // disp(x, "x_2"); 


                            g = g - alpha_f*Ap;
                            this->get_fi(x, g, *lb, *ub, fi);

                            help_f1 = x - (alpha_bar * fi);
                            this->get_projection(help_f1, *lb, *ub, x);

                            A.apply(x, Ax);
                            g = Ax - rhs;
                            this->get_fi(x, g, *lb, *ub, p);
                        }
                    }
                    else
                    {
                        std::cout<<"-------- NO beta beta \n"; 
                        A.apply(beta, Abeta);
                        alpha_cg = dot(g, beta)/dot(beta, Abeta);

                        std::cout<<"delta: "<< delta << "  \n"; 
                        alpha_cg = device::min(alpha_cg, Scalar(delta/norm_infty(beta)));

                        // std::cout<<"alpha_cg: "<< alpha_cg << "  \n"; 
                        // disp(beta,  "beta"); 


                        x = x - alpha_cg*beta;
                        // disp(x, "x_3"); 
                        g = g - alpha_cg*Abeta;

                        this->get_fi(x, g, *lb, *ub, p);
                        // disp(x, "x_4"); 
                    }

                    this->get_fi(x, g, *lb, *ub, fi);
                    this->get_beta(x, g, *lb, *ub, beta);

                    gp = fi+beta;

                    dots(beta, beta, beta_beta, fi, fi, fi_fi, gp, gp, gnorm);

                    gnorm = std::sqrt(gnorm);
                    it++;

                    if(this->verbose()){
                        PrintInfo::print_iter_status(it, {gnorm});
                    }

                    converged = this->check_convergence(it, gnorm, 1, 1);
                }


                // std::cout<<"||x||: "<< norm_infty(x) << "  \n"; 
                // std::cout<<"||ub||: "<< norm_infty(*ub) << "  \n"; 
                // std::cout<<"||lb||: "<< norm_infty(*lb) << "  \n"; 

                // disp(x);
                // exit(0);

                // TODO:: investigate whats wrong with
                this->make_iterate_feasible(x);

                // //cudaProfilerStop();
                // UTOPIA_NO_ALLOC_END();
                return true;
            }


        public:

            void get_fi(const Vector &x, const Vector &g, const Vector &lb, const Vector &ub, Vector & fi) const
            {
                assert(!empty(fi));

                {
                    auto d_lb = const_device_view(lb);
                    auto d_ub = const_device_view(ub);
                    auto d_x  = const_device_view(x);
                    auto d_g  = const_device_view(g);

                    parallel_each_write(fi, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        Scalar li = d_lb.get(i);
                        Scalar ui = d_ub.get(i);
                        Scalar xi = d_x.get(i);
                        Scalar gi = d_g.get(i);

                        if(li < xi && xi < ui){
                            return gi;
                        }
                        else{
                            return 0.0;
                        }

                    });

                }
            }

            Scalar get_alpha_f(const Vector &x, const Vector &p, const Vector &lb, const Vector &ub, Vector & help_f1, Vector & help_f2) const
            {
                assert(!empty(help_f1));
                assert(!empty(help_f2));

                {
                    auto d_lb = const_device_view(lb);
                    auto d_ub = const_device_view(ub);
                    auto d_x  = const_device_view(x);
                    auto d_p  = const_device_view(p);

                    parallel_each_write(help_f1, UTOPIA_LAMBDA(const SizeType i) -> Scalar
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

                    parallel_each_write(help_f2, UTOPIA_LAMBDA(const SizeType i) -> Scalar
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

                // Vector helpp1 = x-lb;
                // Vector helpp2 = ub-x;
                // Vector mu = e_mul(helpp1, helpp2); 

                // Scalar mu0 = min(mu); 
                // std::cout<<"mu0: "<< mu0 << " \n";


                // return std::min(mu0, multi_min(help_f1, help_f2));                

                return multi_min(help_f1, help_f2);
            }


            void get_beta(const Vector &x, const Vector &g, const Vector &lb, const Vector &ub, Vector & beta) const
            {
                assert(!empty(beta));

                {
                    auto d_lb = const_device_view(lb);
                    auto d_ub = const_device_view(ub);
                    auto d_x  = const_device_view(x);
                    auto d_g  = const_device_view(g);

                    parallel_each_write(beta, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                    {
                        Scalar li = d_lb.get(i);
                        Scalar ui = d_ub.get(i);
                        Scalar xi = d_x.get(i);
                        Scalar gi = d_g.get(i);

                        if(device::abs(li -  xi) < 1e-14)
                        {
                            return device::min(0.0, gi);
                        }
                        else if(device::abs(ui -  xi) < 1e-14)
                        {
                            return device::max(0.0, gi);
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
                // Vector y_old;
                // Vector y = local_values(n_loc, 1.0);
                assert(!empty(help_f2));

                // if(empty(help_f2))
                //     help_f2 = local_values(n_loc, 1.0);
                // else
                help_f2.set(1.0);

                SizeType it = 0;
                bool converged = false;
                Scalar gnorm, lambda = 0.0, lambda_old;

                while(!converged)
                {
                    help_f1 = help_f2;
                    A.apply(help_f1, help_f2);
                    help_f2  = Scalar(1.0/Scalar(norm2(help_f2)))*help_f2;

                    lambda_old = lambda;

                    A.apply(help_f2, help_f1);
                    lambda = dot(help_f2, help_f1);

                    fi = help_f2 - help_f1;
                    gnorm = norm2(fi);

                    converged  = ((gnorm < eps_eig_est_) || (std::abs(lambda_old-lambda) < eps_eig_est_) || it > power_method_max_it_) ?  true: false;

                    it=it+1;
                }

                if(this->verbose())
                    std::cout<<"Power method converged in "<< it << " iterations. Largest eig: "<< lambda << "  \n";

                return lambda;
            }

        public:
            void init_memory(const SizeType & ls) override
            {
                auto zero_expr = local_zeros(ls);

                fi = zero_expr;
                beta = zero_expr;
                gp = zero_expr;
                p = zero_expr;
                y = zero_expr;
                Ap = zero_expr;
                Abeta = zero_expr;
                Ax = zero_expr;
                g = zero_expr;
                help_f1 = zero_expr;
                help_f2 = zero_expr;

                initialized_ = true;
                loc_size_ = ls;
            }


        private:
            Vector fi, beta, gp, p, y, Ap, Abeta, Ax, g, help_f1, help_f2;

            Scalar eps_eig_est_;
            SizeType power_method_max_it_;

            bool initialized_;
            SizeType loc_size_;

    };
}

#endif //UTOPIA_CONSTRAINT_QP_MPRGP
