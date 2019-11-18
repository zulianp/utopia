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

    template<class Matrix, class Vector>
    class MPGRP final:  public MatrixFreeQPSolver<Vector>, public QPSolver<Matrix, Vector>
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

                SizeType loc_size_rhs = local_size(rhs);

                //If it is the first time we are here does not make sense to check the sizes
                if(!initialized_ || !rhs.comm().conjunction(loc_size_ == loc_size_rhs)) {
                    init(loc_size_rhs);
                }

                return aux_solve(A, rhs, sol, box);
            }

            bool solve(const Matrix &A, const Vector &rhs, Vector &sol) override
            {
                // auto A_op_ptr = utopia::op_ref(A);
                this->fill_empty_bounds();
                auto &box = this->get_box_constraints();

                SizeType loc_size_rhs = local_size(rhs);

                //If it is the first time we are here does not make sense to check the sizes
                if(!initialized_ || !rhs.comm().conjunction(loc_size_ == loc_size_rhs)) {
                    init(loc_size_rhs);
                }

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

                A.apply(x, Ax);
                g = Ax - rhs;

                this->get_fi(x, g, *lb, *ub, fi);
                this->get_beta(x, g, *lb, *ub, beta);

                gp = fi + beta;
                p = fi;

                dots(beta, beta, beta_beta, fi, fi, fi_fi);

                while(!converged)
                {
                    if(beta_beta <= (gamma*gamma * fi_fi))
                    {
                        A.apply(p, Ap);

                        dots(p, Ap, pAp, g, p, gp_dot);

                        alpha_cg = gp_dot/pAp;
                        y = x - alpha_cg*p;
                        alpha_f = get_alpha_f(x, p, *lb, *ub, help_f1, help_f2);

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

                            help_f1 = x - (alpha_bar * fi);
                            this->get_projection(help_f1, *lb, *ub, x);

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

                    dots(beta, beta, beta_beta, fi, fi, fi_fi, gp, gp, gnorm);

                    gnorm = std::sqrt(gnorm);
                    it++;

                    if(this->verbose()){
                        PrintInfo::print_iter_status(it, {gnorm});
                    }

                    converged = this->check_convergence(it, gnorm, 1, 1);
                }

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


            void init(const SizeType &ls)
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
