#ifndef UTOPIA_CONSTRAINT_QP_MPRGP
#define UTOPIA_CONSTRAINT_QP_MPRGP

#include <string>
#include "utopia_BoxConstraints.hpp"
#include "utopia_QPSolver.hpp"

namespace  utopia
{

    template<class Matrix, class Vector>
    class MPGRP final:  public MatrixFreeQPSolver<Vector>, public QPSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                   Scalar;
        typedef utopia::LinearSolver<Matrix, Vector>    Solver;

        public:
            MPGRP()
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


                // init(local_size(rhs).get(0));
                return aux_solve(A, rhs, sol, box);
            }

            bool solve(const Matrix &A, const Vector &rhs, Vector &sol) override
            {
                auto A_op_ptr = utopia::op_ref(A);
                this->fill_empty_bounds(); 
                auto &box = this->get_box_constraints();

                // init(local_size(rhs).get(0));
                return aux_solve(*A_op_ptr, rhs, sol, box);
            }


        private:
            bool aux_solve(const Operator<Vector> &A, const Vector &rhs, Vector &x, const BoxConstraints<Vector> & constraints)
            {
                const auto &ub = constraints.upper_bound();
                const auto &lb = constraints.lower_bound();

                if(this->verbose()){
                    this->init_solver("MPGRP", {"it", "|| g ||"});
                }

                const Scalar gamma =1; 
                const Scalar alpha_bar = 1.95/this->get_normA(); 

                SizeType it =0; 
                bool converged= false; 
                Scalar gnorm; 

                Scalar alpha_cg, alpha_f, beta_sc; 

                Vector Px; 
                this->get_projection(x, *lb, *ub, Px); 
                x = Px; 

                Vector Ax;
                A.apply(x, Ax);
                Vector g = Ax - rhs; 

                Vector fi, beta, gp, p, y, Ap, Abeta; 
                this->get_fi(x, g, *lb, *ub, fi); 
                this->get_beta(x, g, *lb, *ub, beta); 
                gp = fi + beta; 
                p = fi; 


                while(!converged)
                {
                    if(dot(beta, beta) <= (gamma*gamma * dot(fi,fi)))
                    {
                        A.apply(p, Ap);
                        alpha_cg = dot(g, p)/dot(p,Ap);
                        y = x - alpha_cg*p;
                        alpha_f = get_alpha_f(x, p, *lb, *ub);

                        if(alpha_cg <= alpha_f)
                        {
                            x = y; 
                            g = g - alpha_cg*Ap;
                            this->get_fi(x, g, *lb, *ub, fi); 
                            beta_sc = dot(fi,Ap)/dot(p,Ap);
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

                    it++; 
                    
                    gnorm = norm2(gp); 
                    converged = (gnorm < this->atol() || it > this->max_it()) ? true : false; 

                    if(this->verbose())
                        PrintInfo::print_iter_status(it, {gnorm});

                }

    
                return true;
            }

          

        private:

            void get_fi(const Vector &x, const Vector &g, const Vector &lb, const Vector &ub, Vector & fi) const
            {
                fi = local_values(local_size(x).get(0), 0); 

                {
                    Read<Vector> r_ub(ub), r_lb(lb), r_x(x), r_g(g);
                    
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
                Vector alpha_f1 = local_values(local_size(x).get(0), 1e15); 
                Vector alpha_f2 = local_values(local_size(x).get(0), 1e15); 

                {
                    Read<Vector> r_ub(ub), r_lb(lb), r_x(x), r_g(p);
                    
                    each_write(alpha_f1, [&lb, &x, &p](const SizeType i) -> double 
                    {
                        Scalar li = lb.get(i); 
                        Scalar xi = x.get(i); 
                        Scalar pi = p.get(i);

                        if(pi > 0)
                        {
                            return (xi-li)/pi; 
                        }
                        else
                        {
                            return 1e15; 
                        }

                    }   );



                    each_write(alpha_f2, [&ub, &x, &p](const SizeType i) -> double 
                    {
                        Scalar ui = ub.get(i); 
                        Scalar xi = x.get(i); 
                        Scalar pi = p.get(i);

                        if(pi < 0)
                        {
                            return (xi-ui)/pi; 
                        }
                        else
                        {
                            return 1e15; 
                        }

                    }   );

                }

                Scalar alpha_f2_min = min(alpha_f2); 
                Scalar alpha_f1_min = min(alpha_f1); 

                return std::min(alpha_f1_min, alpha_f2_min); 
            }


            void get_beta(const Vector &x, const Vector &g, const Vector &lb, const Vector &ub, Vector & beta) const
            {
                beta = local_values(local_size(x).get(0), 0.0); 

                {
                    Read<Vector> r_ub(ub), r_lb(lb), r_x(x), r_g(g);
                    
                    each_write(beta, [&ub, &lb, &x, &g](const SizeType i) -> double 
                    {
                        Scalar li = lb.get(i); 
                        Scalar ui = ub.get(i); 
                        Scalar xi = x.get(i); 
                        Scalar gi = g.get(i);

                        if(std::abs(li -  xi) < 1e-13)
                        {
                            return std::min(0.0, gi); 
                        }
                        else if(std::abs(ui -  xi) < 1e-13)
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


            Scalar get_normA()
            {
                // return 36; 
                return 1000; 
            }


        // private:
        //     Vector r, q, d, Hd; // Speed of algo, can be improve big time... 

    };
}

#endif //UTOPIA_CONSTRAINT_QP_MPRGP
