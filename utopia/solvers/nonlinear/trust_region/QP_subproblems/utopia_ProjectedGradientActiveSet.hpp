#ifndef TR_PROJECTED_GRAD_ACTIVE_SET_SUBPROBLEM
#define TR_PROJECTED_GRAD_ACTIVE_SET_SUBPROBLEM
#include <string>
#include "utopia_TRSubproblem.hpp"
#include "utopia_BoxConstraints.hpp"

namespace  utopia 
{

    template<class Matrix, class Vector>
    class ProjectedGradientActiveSet : public TRBoxSubproblem<Matrix, Vector>, 
                                   public MatrixFreeSolverInterface<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;

        typedef utopia::LinearSolver<Matrix, Vector>            LinearSolver;
        typedef utopia::TRBoxSubproblem<Matrix, Vector>         TRBoxSubproblem;


        public:
            ProjectedGradientActiveSet(  const Parameters params = Parameters()):
                                        TRBoxSubproblem(params)
            {
                
            }
            
            virtual ~ProjectedGradientActiveSet( ){}

            ProjectedGradientActiveSet * clone() const override
            {
                return new ProjectedGradientActiveSet(*this);
            }


            virtual bool tr_constrained_solve(const Matrix &H, const Vector &g, Vector &s, const BoxConstraints<Vector> & constraints) override
            {
                return inner_solve(g, s, constraints, H); 
            };

            virtual bool tr_constrained_solve(const Vector &g, Vector &s, const BoxConstraints<Vector> & constraints) override
            {
                return inner_solve(g, s, constraints); 
            }


        private:


            bool inner_solve(const Vector &g, Vector &s, const BoxConstraints<Vector> & constraints, const Matrix &H = Matrix() )
            {
                bool approx_flg = false; 
                if(empty(H))
                    approx_flg = true; 


                const auto &ub = constraints.upper_bound();
                const auto &lb = constraints.lower_bound();


                cp_.tr_constrained_solve(H, g, s, constraints); 

                Vector feasible_set; 
                
                Vector help_g = H*s;    // this->apply_H(s, help_g); 
                Vector grad_qp_fun = g + help_g;  // minus is missing 

                // building feasible set 
                this->build_feasible_set(s, *ub, *lb, feasible_set); 
                SizeType feasible_variables = sum(feasible_set); 


            
                if(feasible_variables == 0) // all variables are feasible => perform Newton step on whole matrix
                {
                    std::cout<<"------- Feasible set empty => return CP point --------- \n";     
                    return false; 
                }
                else if(size(feasible_set).get(0)==feasible_variables)
                {
                    std::cout<<"------- Active set empty => apply inverse --------- \n";     
                    Vector  local_corr; 
                    auto linear_solver = std::make_shared<Factorization<Matrix, Vector> >(); // this->apply_Hinv(grad_quad_fun, local_corr); 
                    linear_solver->solve(H, -1.0 * grad_qp_fun, local_corr); 

                    Vector lb_new = *lb - s; 
                    Vector ub_new = *ub - s; 

                    // alternative solution is to do full Newton step and then make iterate feasible by projecting back to the boundary
                    Scalar alpha_star = compute_alpha_max(0*s, lb_new, ub_new, local_corr); 
                    // std::cout<<"alpha_star: "<< alpha_star << " \n"; 
                    s += alpha_star * local_corr;    
                    
                }
                else
                {
                    // std::cout<<"------- Projected CG method --------- \n"; 
                    Vector  local_corr = 0*s; 

                    Vector lb_new = *lb - s; 
                    Vector ub_new = *ub - s; 
                    projected_cg(grad_qp_fun, feasible_set, ub_new, lb_new, H, local_corr); 

                    s += local_corr; 
                }

                return true; 
            }



            void projected_cg(const Vector & grad, const Vector & feasible_set, const Vector &ub, const Vector &lb, const Matrix & H, Vector & x)
            {
                Scalar r_norm0 = norm2(grad); 
                Scalar alpha, betta, alpha_max; 

                Vector r = grad; 
                Vector g;
                get_projection(r, feasible_set,  g); 

                Vector d = -1.0 * g; 
                Vector r_plus, g_plus; 

                Scalar d_norm = norm2(d); 

                // disp(feasible_set); 

                // std::cout<<"d_norm: "<< d_norm << "  \n"; 

                if(d_norm< 1e-15 || !std::isfinite(d_norm))
                    return; 

                SizeType max_it = 200, it = 0; 
                bool converged = false; 

                while(!converged)
                {
                    alpha = dot(r,g)/dot(d, H*d); 

                    if(!std::isfinite(alpha))
                        return; 

                    alpha_max = compute_alpha_max(x, lb, ub, d); 

                   //  std::cout<<"alpha: "<< alpha << "  alpha_max: "<< alpha_max << "  \n"; 

                    if(alpha > alpha_max)
                    {
                        x += alpha_max * d; 
                        return; 
                    }

                    x +=  alpha *d; 

                    r_plus = r + (alpha * H * d); 
                    get_projection(r_plus, feasible_set, g_plus); 


                    betta = dot(r_plus, g_plus)/dot(r, g); 

                    if(!std::isfinite(betta))
                        return; 

                    d = betta* d - g_plus; 
                    g = g_plus; 
                    r = r_plus; 

                    converged = (it > max_it || norm2(r) < std::min(0.1, std::sqrt(r_norm0)) * r_norm0 || alpha < 1e-12 || betta < 1e-12) ? true : false; 
                    // converged = (it > max_it || norm2(r) < 1e-5 )? true : false; 

                    it++; 

                   // std::cout<<"  it: "<< it <<  "  alpha: "<<alpha <<  "  alpha_max: "<<alpha_max <<"   betta:  "<< betta << " ||r||  "<<  norm2(r) << "  \n"; 

                }

                // exit(0); 
            }



            void get_projection(const Vector & r,  const Vector & feasible_set, Vector & Pr)
            {
                // simple identity preconditioner - orthogonal to the nullspace of the equality constraints 
                Pr = e_mul(r, feasible_set); 
            }



        private:
            SizeType get_local_size_feasible_set(const Vector & feasible_set) const
            {
                SizeType local_feasible_set = 0; 
                each_read(feasible_set, [&local_feasible_set](const SizeType i, const Scalar value) 
                {
                    if(approxeq(value, 1.0)) 
                        local_feasible_set++; 
                });

                return local_feasible_set; 
            }


            void build_feasible_set(const Vector &x, const Vector &ub, const Vector &lb, Vector & feasible_set) const 
            { 
                if(empty(feasible_set) || local_size(feasible_set)!=local_size(x))
                    feasible_set = local_zeros(local_size(x)); 

                Read<Vector> r_ub(ub), r_lb(lb), r_x(x);
                each_write(feasible_set, [ub, lb, x](const SizeType i) -> double { 
                            Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);  
                            if(li < xi && xi < ui)
                                return 1.0; 
                            else
                                return 0.0; }   );
            }


            void build_active_set(const Vector &x, const Vector &ub, const Vector &lb, Vector & active_set) const 
            { 
                if(empty(active_set) || local_size(active_set)!=local_size(x))
                    active_set = local_zeros(local_size(x)); 

                Read<Vector> r_ub(ub), r_lb(lb), r_x(x);
                each_write(active_set, [ub, lb, x](const SizeType i) -> double { 
                            Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);  
                            if(li < xi && xi < ui)
                                return 0.0; 
                            else
                                return 1.0; }   );
            }



        // to be rechecked tomorrow 
        Scalar compute_alpha_max(const Vector & x, const Vector & lb, const Vector & ub, const Vector & d) const
        {
            Vector alpha_stars = local_values(local_size(x).get(0), 1.0);  

            {
                Read<Vector>  rd(d); 
                Read<Vector>  ru(ub); 
                Read<Vector>  rlb(lb); 
                Read<Vector>  rs(x); 

                Write<Vector> wv(alpha_stars); 

                auto rr = range(alpha_stars); 

                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {                
                    Scalar val = 0.0; 
                    Scalar di = d.get(i);  Scalar xi = x.get(i); Scalar li = lb.get(i);  Scalar ui = ub.get(i); 

                    if(di!=0.0)
                        val = (di>0) ? (ui - xi)/di : (li - xi)/di; 

                    // checks for nans and infs, also alpha needs to be positive... 
                    if(val < 1.0 && std::isfinite(val) && val > 0.0)
                        alpha_stars.set(i, val); 
                }     
            }

            return min(alpha_stars); 
        }






        private:  
            GeneralizedCauchyPoint<Matrix, Vector> cp_; 
        
    };
}

#endif //TR_PROJECTED_GRAD_ACTIVE_SET_SUBPROBLEM
