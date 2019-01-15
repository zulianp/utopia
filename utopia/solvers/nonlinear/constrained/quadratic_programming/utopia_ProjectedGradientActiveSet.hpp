#ifndef TR_PROJECTED_GRAD_ACTIVE_SET_SUBPROBLEM
#define TR_PROJECTED_GRAD_ACTIVE_SET_SUBPROBLEM
#include <string>
#include "utopia_BoxConstraints.hpp"
#include "utopia_QPSolver.hpp"

namespace  utopia 
{

    template<class Matrix, class Vector>
    class ProjectedGradientActiveSet final:  public MatrixFreeQPSolver<Vector>, public QPSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                   Scalar;
        typedef utopia::LinearSolver<Matrix, Vector>    Solver;

        public:
            ProjectedGradientActiveSet()
            {
                
            }
            
            void read(Input &in) override
            {
                MatrixFreeQPSolver<Vector>::read(in);
                QPSolver<Matrix, Vector>::read(in);

                in.get("Cauchy-point", cp_);

                if(precond_) {
                    in.get("precond", *precond_);
                }
                if(linear_solver_) {
                    in.get("linear-solver", *linear_solver_);
                }                                                
            }


            void print_usage(std::ostream &os) const override
            {
                MatrixFreeQPSolver<Vector>::print_usage(os);
                QPSolver<Matrix, Vector>::print_usage(os);

                this->print_param_usage(os, "Cauchy-point", "GeneralizedCauchyPoint", "Input parameters for Generalized Cauchy point solver.", "-"); 
                this->print_param_usage(os, "precond", "Preconditioner", "Input parameters for Preconditioner.", "-"); 
                this->print_param_usage(os, "linear-solver", "LinearSolver", "Input parameters for LinearSolver.", "-"); 
            }



            ProjectedGradientActiveSet * clone() const override
            {
                return new ProjectedGradientActiveSet(*this);
            }

            void set_preconditioner(const std::shared_ptr<Preconditioner<Vector> > &precond)
            {
                precond_ = precond;
            }
            
            bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) override 
            {
                auto &box = this->get_box_constraints(); 
                init(local_size(rhs).get(0));
                return aux_solve(A, -1.0*rhs, sol, box); 
            }

            bool solve(const Matrix &A, const Vector &rhs, Vector &sol) override
            {
                auto A_op_ptr = utopia::op_ref(A);
                auto &box = this->get_box_constraints(); 
                return aux_solve(*A_op_ptr, -1.0 *rhs, sol, box); 
            }


            void set_linear_solver(const std::shared_ptr<Solver> &linear_solver)
            {
                linear_solver_ = linear_solver; 
            }

        private:
            bool aux_solve(const Operator<Vector> &H, const Vector &g, Vector &s, const BoxConstraints<Vector> & constraints)
            {
                const auto &ub = constraints.upper_bound();
                const auto &lb = constraints.lower_bound();

                if(this->verbose())
                    this->init_solver("ProjectedGradientActiveSet", {" "});

                cp_.set_box_constraints(constraints); 
                cp_.solve(H, -1.0*g, s); 

                Vector feasible_set; 
                
                Vector help_g;
                H.apply(s, help_g); 

                Vector grad_qp_fun = g + help_g;  // minus is missing 

                // building feasible set 
                this->build_feasible_set(s, *ub, *lb, feasible_set); 
                SizeType feasible_variables = sum(feasible_set); 

            
                if(feasible_variables == 0) // no feasible variables => return CP
                {
                    return false; 
                }
                else if(size(feasible_set).get(0)==feasible_variables) // all variables are feasible, use Newton's method to solve the system
                {
                    Vector  local_corr = local_zeros(local_size(s).get(0)); 
                    if(const MatrixOperator<Matrix, Vector> * H_matrix = dynamic_cast<const MatrixOperator<Matrix, Vector> *>(&H))
                    {
                        if(linear_solver_)
                        {
                            auto linear_solver = std::make_shared<Factorization<Matrix, Vector> >(); 
                            linear_solver->solve(*(H_matrix->get_matrix()), -1.0 * grad_qp_fun, local_corr); 
                        }
                        else if(precond_)
                            precond_->apply(-1.0 * grad_qp_fun, local_corr);                             
                    }
                    else
                    {
                        if(precond_)
                            precond_->apply(-1.0 * grad_qp_fun, local_corr); 
                    }

                    Vector lb_new = *lb - s; 
                    Vector ub_new = *ub - s; 

                    // alternative solution is to do full Newton step and then make iterate feasible by projecting back to the boundary
                    Scalar alpha_star = compute_alpha_max(0*s, lb_new, ub_new, local_corr); 
                    // std::cout<<"alpha_star: "<< alpha_star << " \n"; 
                    s += alpha_star * local_corr;   
                    
                }
                else
                {
                    Vector  local_corr = 0*s; 

                    Vector lb_new = *lb - s; 
                    Vector ub_new = *ub - s; 
                    projected_cg(H, grad_qp_fun, feasible_set, ub_new, lb_new, local_corr); 

                    s += local_corr; 
                }

                return true; 
            }


            void projected_cg(const Operator<Vector> &H, const Vector & grad, const Vector & feasible_set, const Vector &ub, const Vector &lb, Vector & x)
            {
                Scalar r_norm0 = norm2(grad), r_norm; 
                Scalar alpha, betta, alpha_max, rg_old, d_norm; 

                r = grad; 
                get_projection(r, feasible_set,  q); 
                d = -1.0 * q; 

                d_norm = norm2(d); 

                if(d_norm< this->stol()|| !std::isfinite(d_norm))
                    return; 

                if(this->verbose())
                    this->init_solver("Projected CG - inner solver", {"it", "alpha", "alpha_max", "betta", "|| r ||"});

                SizeType it = 0; 
                bool converged = false; 

                while(!converged)
                {
                    H.apply(d, Hd); 

                    rg_old = dot(r, q); 
                    alpha = rg_old/dot(d, Hd); 

                    if(!std::isfinite(alpha))
                        return; 

                    alpha_max = compute_alpha_max(x, lb, ub, d); 
                    if(alpha > alpha_max)
                    {
                        x += alpha_max * d; 
                        return; 
                    }

                    x +=  alpha * d; 
                    r += alpha * Hd;

                    get_projection(r, feasible_set, q); 
                    betta = dot(r, q)/rg_old; 

                    if(!std::isfinite(betta))
                        return; 

                    d = betta* d - q; 

                    r_norm = norm2(r); 

                    // following Conn, Gould, Toint... 
                    converged = (it > this->max_it() || r_norm < std::min(0.1, std::sqrt(r_norm0)) * r_norm0 || alpha < 1e-12 || betta < 1e-12) ? true : false; 
                    // converged = (it > max_it || norm2(r) < 1e-5 )? true : false; 

                    it++; 

                    if(this->verbose())
                        PrintInfo::print_iter_status(it, {alpha, alpha_max, betta, r_norm});
                }

            }



            void get_projection(const Vector & r,  const Vector & feasible_set, Vector & Pr)
            {
                // simple identity preconditioner - orthogonal to the nullspace of the equality constraints 
                Pr = e_mul(r, feasible_set); 
                // TODO:: extend to other possibilities...
            }


            void set_GCP_memory(const SizeType & m)
            {
                cp_.set_memory_size(m); 
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

                {
                    Read<Vector> r_ub(ub), r_lb(lb), r_x(x);
                    each_write(feasible_set, [&ub, &lb, &x](const SizeType i) -> double { 
                                Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);  
                                if(li < xi && xi < ui)
                                    return 1.0; 
                                else
                                    return 0.0; }   );
                }
            }


            void build_active_set(const Vector &x, const Vector &ub, const Vector &lb, Vector & active_set) const 
            { 
                if(empty(active_set) || local_size(active_set)!=local_size(x))
                    active_set = local_zeros(local_size(x)); 
                {
                    Read<Vector> r_ub(ub), r_lb(lb), r_x(x);
                    each_write(active_set, [&ub, &lb, &x](const SizeType i) -> double { 
                                Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);  
                                if(li < xi && xi < ui)
                                    return 0.0; 
                                else
                                    return 1.0; }   );
                }
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

            void init(const SizeType &ls)
            {
                auto zero_expr = local_zeros(ls);

                //resets all buffers in case the size has changed
                if(!empty(r)) {
                    r = zero_expr;
                }

                if(!empty(q)) {
                    q = zero_expr;
                }

                if(!empty(d)) {
                    d = zero_expr;
                }

                if(!empty(Hd)) {
                    Hd = zero_expr;
                }
            }


        private:  
            GeneralizedCauchyPoint<Matrix, Vector> cp_; 
            std::shared_ptr<Preconditioner<Vector> > precond_;
            std::shared_ptr<Solver> linear_solver_;    

            Vector r, q, d, Hd; 
        
    };
}

#endif //TR_PROJECTED_GRAD_ACTIVE_SET_SUBPROBLEM
