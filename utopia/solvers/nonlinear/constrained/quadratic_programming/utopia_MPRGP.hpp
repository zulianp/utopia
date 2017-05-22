/*
* @Author: Alena Kopanicakova
* @Date:   2017-05-19
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-05-22
*/

#ifndef UTOPIA_QPSOLVER_BOX_MPRGP_NEWTON_HPP
#define UTOPIA_QPSOLVER_BOX_MPRGP_NEWTON_HPP

#include "utopia_Wrapper.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_BoxConstraints.hpp"
#include <vector>
#include <limits>



namespace utopia 
{
    template<class Matrix, class Vector>
    class MPRGP : public IterativeSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)    SizeType;

        typedef utopia::BoxConstraints<Vector> BoxConstraints;
        
    public:
        
        MPRGP(  const Parameters params = Parameters())
        {
            set_parameters(params);
        }


        bool solve(const Matrix &A, const Vector &b, Vector &x)  override
        {
           using namespace utopia;

           SizeType it = 0; 
           bool converged = false;

           Scalar normA, alpha_cg, alpha_f, g_norm = 9999, beta2; 
           gersgorin(A, normA); 

           Scalar gamma = 1.0, alpha_bar = 1.95/normA; 


           Vector Px =x; 
           Vector lb = *constraints_->lower_bound(); 
           Vector ub = *constraints_->upper_bound(); 

           get_projection(x, lb, ub, Px); 

           // TODO:: check this one 
           x = Px; 

           Vector g = A * x - b; 
           Vector fi, beta; 

           get_fi(x, g, lb, ub, fi); 
           get_beta(x,  g, lb, ub, beta); 

           
           Vector gp = fi + beta; 
           Vector p = fi; 

           Vector Ap, y; 

            if(this->verbose())
            {
                this->init_solver("MPRGP", {" it. ", "|| g ||"});
                PrintInfo::print_iter_status(it, {g_norm});
            }
            
            it++;

            while(!converged)
            {

                // 1. Proportional x_k. Trial conjugate gradient step
                if(dot(beta,beta) <= gamma * gamma * dot(fi,fi))
                {
                    Ap = A*p; 
                    alpha_cg = dot(g,p)/dot(p,Ap); 
                    y = x - alpha_cg * p; 
                    get_alpha_f(x, p, lb, ub, alpha_f); 

                    if(std::min(alpha_cg, alpha_f) == std::numeric_limits<double>::infinity())
                    {
                        std::cout<<"MPRGP:: no solution ! \n"; 
                        return 0; 
                    }

                    // 2. Conjugate gradient step
                    if(alpha_cg <= alpha_f)
                    {
                        x = y; 
                        g = g - alpha_cg * Ap; 
                        get_fi(x, g, lb, ub, fi); 
                        beta2 = dot(fi, Ap)/dot(p, Ap); 
                        p = fi - beta2 * p; 
                    }
                    // 3. Expansion step
                    else
                    {  
                        x = x - alpha_f * p; 
                        g = g - alpha_f * Ap; 
                        get_fi(x, g, lb, ub, fi); 
                        x = x - alpha_bar*fi;
                        get_projection(x, lb, ub, Px); 
                        x = Px; 
                        g = A*x - b; 
                        get_fi(x, g, lb, ub, p); 
                    }
                }

                else
                {    
                    // 4. Proportioning step
                    Vector Abeta =  A * beta; 
                    alpha_cg = dot(g, beta)/ dot(beta, Abeta); 
                    x = x - alpha_cg * beta;
                    g = g - alpha_cg * Abeta; 
                    get_fi(x, g, lb, ub, p); 
                }


                get_fi(x, g, lb, ub, fi); 
                get_beta(x,  g, lb, ub, beta); 

                gp = fi + beta; 

                g_norm = norm2(gp); 


                // print iteration status on every iteration
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {g_norm});

                // // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 1, 1);
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
        std::shared_ptr<BoxConstraints> constraints_; 

        bool gersgorin(const Matrix & A, Scalar & normA)
        {
            Vector s = diag(A); 

            Matrix D = local_sparse(local_size(A).get(0), local_size(A).get(1), 1); 
            D = diag(s); 

            Matrix R = A - D; 
            Vector v = s; 

            // TODO: we should have abs of entry, not of result... 
            v = sum(R, 1); 
            v = abs(v); 

            v += s; 

            normA = max(v); 

            return true; 
        }


        // * @brief      Function projects constraints, such that \f$ lb < x < ub \f$
        bool get_projection(const Vector & x, const Vector &lb, const Vector &ub, Vector & Px)
        {
            Px = x; 

            {
                Read<Vector> r_ub(ub), r_lb(lb);
                each_transform(x, Px, [&ub, &lb](const SizeType i, const Scalar entry) -> double  
                { 
                    Scalar ui = ub.get(i), li = lb.get(i); 
                    if(entry > ui && entry < li)
                    {
                        return (std::max(li, std::min(ui, entry)));
                    }
                    else
                        return entry;   
                }    );
            }
            return true;
        }




        bool get_fi(const Vector & x, const Vector & g, const Vector &lb, const Vector &ub, Vector & fi)
        {
            fi = local_zeros(local_size(x)); 
            {
                Write<Vector> w (fi);
                Read<Vector> rg(g), rx(x), rl(lb), rb(ub);
                Range rhs_range = range(fi);;
                for (SizeType i = rhs_range.begin(); i != rhs_range.end() ; i++)
                {
                    Scalar ui = ub.get(i), li = lb.get(i), gi = g.get(i), xi = x.get(i); 
                    if(xi > li && xi < ui)
                        fi.set(i, gi);
                }
            }
            return true; 
        }



        bool get_beta(const Vector & x, const Vector & g, const Vector &lb, const Vector &ub, Vector & beta)
        {
            beta = local_zeros(local_size(x)); 
            {
                Write<Vector> w (beta);
                Read<Vector> rg(g), rx(x), rl(lb), rb(ub);
                Range beta_range = range(beta);;
                for (SizeType i = beta_range.begin(); i != beta_range.end() ; i++)
                {
                    Scalar ui = ub.get(i), li = lb.get(i), gi = g.get(i), xi = x.get(i); 
                    if(xi==li)
                        beta.set(i, std::min(gi, 0.0) );

                    if(xi==ui)
                        beta.set(i, std::max(gi, 0.0) );
                }
            }
            return true; 
        }




        bool get_alpha_f(const Vector & x, const Vector & p, const Vector &lb, const Vector &ub, Scalar & alpha)
        {

            Scalar alpha1, alpha2;
            Vector a1 = x, a2 = x; 

            {
                Write<Vector> w (a1);
                Read<Vector> rg(p), rx(x), rl(lb);
                Range beta_range = range(a1);;
                for (SizeType i = beta_range.begin(); i != beta_range.end() ; i++)
                {
                    Scalar li = lb.get(i), pi = p.get(i), xi = x.get(i); 
                    if(pi > 0)
                        a1.set(i, (xi - li)/ pi );
                    else
                        a1.set(i, std::numeric_limits<double>::infinity());
                }
            }

            {
                Write<Vector> w (a2);
                Read<Vector> rg(p), rx(x), rl(ub);
                Range beta_range = range(a1);;
                for (SizeType i = beta_range.begin(); i != beta_range.end() ; i++)
                {
                    Scalar ui = ub.get(i), pi = p.get(i), xi = x.get(i); 
                    if(pi < 0)
                        a2.set(i, (xi - ui)/ pi );
                    else
                        a2.set(i, std::numeric_limits<double>::infinity());
                }
            }
            

            alpha1 = min(a1); 
            alpha2 = min(a2); 
            alpha = std::min(alpha1, alpha2); 

            return true; 
        }


    };

}
#endif //UTOPIA_QPSOLVER_BOX_MPRGP_NEWTON_HPP
