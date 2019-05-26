#ifndef UTOPIA_ASTRUM_HPP
#define UTOPIA_ASTRUM_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{
    // since $F = - \nabla f(x)$
    template<class Matrix, class Vector>
    class DAEFormFunction final: public Function<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::Function<Matrix, Vector>    Fun;

        DAEFormFunction(const std::shared_ptr<Function<Matrix, Vector> > & fun):
        fun_(fun)
        {

        }

        bool value(const Vector &/*point*/, Scalar & value) const override
        {
            value = 0; 
            return false; // should not be necessary
        }

        bool gradient(const Vector & x, Vector &g) const override
        {
            fun_->gradient(x, g); 
            g = -1.0 * g; 
            return true; 
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            fun_->hessian(x, H); 
            H = -1.0 * H; 
            return true; 
        }

        private:
            std::shared_ptr<Fun> fun_;   
    };


    template<class Matrix, class Vector>
    class ASTRUM : public NewtonBase<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        typedef typename NewtonBase<Matrix, Vector>::Solver Solver;
        typedef utopia::LSStrategy<Vector>                  LSStrategy; 

        using NewtonBase<Matrix, Vector>::print_statistics; 


    public:
       ASTRUM(  const std::shared_ptr <Solver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >()):
                NewtonBase<Matrix, Vector>(linear_solver), 
                tau_max_(1e14),
                tau_min_(1e-9), 
                tau_zero_user_(-1),
                alpha_treshold_(1e-10), 
                max_inner_it_(5), 
                reset_mass_matrix_(false), 
                is_identity_(false), 
                scaling_user_provided_(false), 
                scaling_(true)
                {
                    verbosity_level_ =  VERBOSITY_LEVEL_NORMAL; 
                }


        void read(Input &in) override
        {
            NewtonBase<Matrix, Vector>::read(in);
            in.get("tau_max", tau_max_);
            in.get("tau_min", tau_min_);
            in.get("tau_zero_user", tau_zero_user_);
            in.get("alpha_treshold", alpha_treshold_);
            in.get("max_inner_it", max_inner_it_);
            in.get("reset_mass_matrix", reset_mass_matrix_); 
            in.get("scaling", scaling_); 

        }

        void print_usage(std::ostream &os) const override
        {
            NewtonBase<Matrix, Vector>::print_usage(os); 
            
        }

        void reset_mass_matrix(const bool reset_mass)
        {
            reset_mass_matrix_ = reset_mass; 
        }


        bool solve(Function<Matrix, Vector> &fun_grad, Vector &x) override
        {
           using namespace utopia;

            std::shared_ptr<Function<Matrix, Vector> > fun_grad_ptr_(&fun_grad, [](Function<Matrix, Vector>*){});
            DAEFormFunction<Matrix, Vector> fun(fun_grad_ptr_); 

            Scalar g_norm=0.0, g_norm_old=0.0, s_norm=0.0, L, tau_old, rho; 
            SizeType it_inner = 0; 

            bool AS_form_used; 

            Vector g, s, g_new, r; 
            s = 0 * x; 
            g_new = 0*x; 
            Matrix H; 


            fun.gradient(x, g); 
            g_norm = norm2(g); 

            fun.hessian(x, H); 


            if(empty(I_) || reset_mass_matrix_==true)
            {
                
                // if(this->verbose())
                // {
                //     std::cout<<"mass matrix not set, using Identity matrix ... \n";
                // }

                I_  = local_identity(local_size(H).get(0), local_size(H).get(1)); 
                is_identity_ = true; 
            }

            if(!scaling_)
            {
                D_ = local_identity(local_size(H).get(0), local_size(H).get(1)); 
                D_inv_ = D_; 

                // if(mpi_world_rank()==0){
                //     std::cout<<" init of scaling matrix.... \n"; 
                // }
            }

            if(empty(D_))
            {

                this->update_scaling_matrices(x, x);

                // initialize D to identity 
                // D_ = local_identity(local_size(H).get(0), local_size(H).get(1)); 
                // D_inv_ = D_; 

                // if(mpi_world_rank()==0){
                //     std::cout<<" init of scaling matrix from x .... \n"; 
                // }
            }            


            Scalar tau = (tau_zero_user_ > 0)? tau_zero_user_ : std::max(1., 1./g_norm);


            bool converged = false; 
            SizeType it = 0; 

            if(verbosity_level_ >= VERBOSITY_LEVEL_NORMAL)
            {
                this->init_solver("ASTRUM", {" it. ", "|| F ||", "|| Delta x || ", "tau", "mu",  "AF_form", "inner it"});
                PrintInfo::print_iter_status(it, {g_norm, 0, tau});
            }

            it++; 

            while(!converged)
            {   
                // old version of Peter 
                // Matrix A = I_ - (D_inv_ * (tau * H) * D_); 
                Matrix A = D_inv_  * (I_ -  (tau * H)) * D_; 

                s = 0*x; 
                Vector g_help = D_inv_ * g; 
                this->linear_solve(A, g_help, s);                


                // linear system 
                // r = g - A*s; 


                Vector x_trial = x + (tau * D_* s); 
                s_norm = norm2(D_*s); 

                fun.gradient(x_trial, g_new); 
                g_norm_old = g_norm; 
                g_norm = norm2(g_new); 


                // TBD:: remove for efficiency... 
                Scalar mu = estimate_mu(g_new, g, s, tau, s_norm); 
                // std::cout<<"mu: "<< mu << "\n"; 

                // if((mu > 0) && (s_norm < g_norm_old))
                // {
                //     std::cout<<"inconsistencies in condition for checking the violation of condition... \n"; 
                // }

                // seems that only condition on norm of correction is not sufficient 
                if(s_norm < g_norm_old && it_inner < max_inner_it_)
                // if((mu < 0) && it_inner < max_inner_it_)
                {
                    // check convergence of fixed point iteration 
                    tau_old = tau; 
                    tau = estimate_tau(g_new, g, s, tau, s_norm); 

                    AS_form_used = true;
                    rho = (residual_monotonicity_test(g_new, g) || (std::abs(tau_old - tau) < 1e-1)) ? 1.0 : 0.0; 

                }
                // mu is positive, so lets take PTC update formula
                else
                {
                    // std::cout<<"mu  is positive => 1st Lipschitz condition violated  .... \n"; 

                    // Kelley, Keyes formula
                    // tau = tau * g_norm_old/g_norm; 

                    // original formula, no residual decrease required
                    tau = estimate_tau(g_new, g, s, tau, s_norm); 

                    // Peters suggestion
                    // tau = 0.5* tau * g_norm_old/g_norm; 

                    AS_form_used = false;
                    rho = 1.0; 
                }


                if(rho>0)
                {
                    // std::cout<<"step taken.... \n"; 
                    if(!scaling_user_provided_ && scaling_)
                    {
                        this->update_scaling_matrices(x, x_trial);
                    }

                    x = x_trial; 
                    g = g_new; 
                    fun.hessian(x, H); 
                    it_inner=0; 
                }
                else
                {
                    // std::cout<<"step not taken.... \n"; 
                    g_norm = g_norm_old; 
                    it_inner++; 
                }


                // print iteration status on every iteration
                if(verbosity_level_ >= VERBOSITY_LEVEL_NORMAL){
                    PrintInfo::print_iter_status(it, {g_norm, s_norm, tau, mu, Scalar(AS_form_used), Scalar(it_inner)});
                }


                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 9e9, s_norm);
                it++;
            }

            return true;
        }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        void set_max_inner_it(const SizeType & max_it)
        {
            max_inner_it_ = max_it; 
        }

        VerbosityLevel verbosity_level() const 
        {
            return verbosity_level_; 
        }

        void verbosity_level(const VerbosityLevel & verbose_level )
        {
            verbosity_level_ = this->verbose() ? verbose_level : VERBOSITY_LEVEL_QUIET;  
        }


        Scalar tau_max() const { return tau_max_; }
        Scalar tau_min() const { return tau_min_; }
        Scalar tau_init() const { return tau_zero_user_; }

        void tau_max(const Scalar & tau_max)    {   tau_max_ = tau_max; }
        void tau_min(const Scalar & tau_min)    {   tau_min_ = tau_min; }
        void tau_init(const Scalar & tau_init)  {   tau_zero_user_ = tau_init; }


        void scaling(const bool & flg)
        {
            scaling_ = flg; 
        }

        bool scaling() const
        {
            return scaling_; 
        }


        void set_scaling_matrix(const Matrix & D)
        {
            D_ = D; 

            Vector d = diag(D_); 
            D_inv_ = diag(1.0/d); 

            scaling_user_provided_ = true; 
            scaling_ = true; 
        }                

        void set_mass_matrix(const Matrix & M)
        {
            I_ = M; 

            Matrix Diff  = local_identity(local_size(M).get(0), local_size(M).get(1)); 
            Diff = Diff - I_; 

            if(max(Diff) > 0)
            {
                is_identity_ = false; 
            }
            else
            {
                is_identity_ = true; 
            }
        }


    protected:
        virtual void print_statistics(  const SizeType & it, const Scalar & g_norm, 
                                        const Scalar & tau,  const SizeType & it_inner) 
        {
            auto rmtr_data_path = Utopia::instance().get("af_data_path");
            if(!rmtr_data_path.empty())
            {
                CSVWriter writer; 
                if (mpi_world_rank() == 0)
                {
                    if(!writer.file_exists(rmtr_data_path))
                    {
                        writer.open_file(rmtr_data_path); 
                        writer.write_table_row<std::string>({"it", "g", "tau", "it_inner"}); 
                    }
                    else
                        writer.open_file(rmtr_data_path); 

                    writer.write_table_row<Scalar>({Scalar(it), g_norm, tau, Scalar(it_inner)}); 
                    writer.close_file(); 
                }
            }
        }

        void update_scaling_matrices(const Vector & x_old, const Vector & x_new)
        {
            Vector x_scaling = local_values(local_size(x_old).get(0), 1.0); 

            {   
                Read<Vector>    r1(x_old), r2(x_new); 
                auto tol = alpha_treshold_; 
                each_write(x_scaling, [&x_old, &x_new, tol](const SizeType i) -> double 
                { 
                    return std::max(std::max(std::abs(x_old.get(i)), std::abs(x_new.get(i))), tol); 
                });

            }

            D_ = diag(1./x_scaling); 
            D_inv_ = diag(x_scaling);         

            // if(mpi_world_rank()==0){
            //     std::cout<<" scaling matrix updated .... \n"; 
            // }
        }




    private: 
        Scalar estimate_tau(const Vector & g_new, const Vector & g, const Vector & s, const Scalar & tau, const Scalar & s_norm)
        {   
            if(is_identity_)
            {
                return estimate_tau_alg(g_new, g, s, tau, s_norm); 
            }
            else
            {
                return estimate_tau_PDE(g_new, g, s, tau, s_norm); 
            }
        }


        Scalar estimate_tau_PDE(const Vector & g_new, const Vector & g, const Vector & s, const Scalar & tau, const Scalar & s_norm)
        {   
            Scalar s_norm2 = norm_l2_2(s);

            Scalar nom = dot(D_inv_*g, s) - s_norm2; 
            Scalar denom = dot(D_inv_*g_new, s) - s_norm2; 

            Scalar tau_new = 0.5*tau * std::abs(nom/denom); 
            // bool flg = this->clamp_tau(tau_new); 
            this->clamp_tau(tau_new); 

            // changing initial guess for fixed point iteration
            if(tau_new==tau_min_)
            {
                tau_new = 1./Scalar(2.*norm_l2_2(g)); 
            }
        
            return tau_new; 
        }


        Scalar estimate_tau_alg(const Vector & g_new, const Vector & g, const Vector & s, const Scalar & tau, const Scalar & s_norm)
        {   

            Scalar s_norm2 = norm2(s);

            Scalar nom = dot(s, g - s); 
            Scalar denom = s_norm2 * norm2(g_new - s); 

            Scalar tau_new = 0.5*tau * std::abs(nom/denom); 
            // bool flg = this->clamp_tau(tau_new); 
            this->clamp_tau(tau_new); 

            // changing initial guess for fixed point iteration
            if(tau_new==tau_min_)
            {
                tau_new = 1./Scalar(2.*norm2(g)); 
            }
        
            return tau_new; 
        }

    

        Scalar estimate_tau_inexact(const Vector & g_new, const Vector & g, const Vector & s, const Vector & r,  const Scalar & tau, const Scalar & s_norm)
        {   
            Scalar s_norm2 = norm_l2_2(s);

            Scalar r_norm = dot(r, I_*s);

            // std::cout<<"r_norm: "<< r_norm << "  \n"; 

            Scalar nom = dot(g, s) - s_norm2 - r_norm; 
            Scalar denom = dot(g_new, s) - s_norm2 - r_norm; 

            Scalar tau_new = 0.5*tau * std::abs(nom/denom); 
            // bool flg = this->clamp_tau(tau_new); 
            this->clamp_tau(tau_new); 
            
            return tau_new; 
        }        


        Scalar estimate_mu(const Vector & g_new, const Vector & g, const Vector & s, const Scalar & tau, const Scalar & s_norm)
        {   
            Scalar s_norm2 = norm_l2_2(s);

            Scalar nom = s_norm2 - dot(g, s); 
            Scalar denom = tau * s_norm2;


            // TODO:: verify 
            // Scalar nom = dot(s, s-g); 
            // Scalar denom = tau * norm2(s)*norm2(s);


           // nom = (dot(s, s - g));
           // denom = norm(s)^2 * tau;
           // mu = nom/denom; 


            return nom/denom; 
        }


        Scalar estimate_mu_inexact(const Vector & g_new, const Vector & g, const Vector & s,  const Vector & r, const Scalar & tau, const Scalar & s_norm)
        {   
            Scalar s_norm2 = norm_l2_2(s);
            Scalar r_norm = dot(r, I_*s);

            Scalar nom = dot(g, s) - s_norm2 - r_norm; 
            Scalar denom = tau * s_norm2;

            return nom/denom; 
        }        




        Scalar norm_l2_2(const Vector & s)
        {   
            return dot(s, I_*s); 
        }

        Scalar norm_l2(const Vector & s)
        {   
            return std::sqrt(norm_l2_2(s)); 
        }


        bool clamp_tau(Scalar & tau)
        {
            if(std::isinf(tau) || tau > tau_max_ )
            {
                tau = tau_max_;        
                return true; 
            }
            else if (std::isnan(tau))
            {
                tau = tau_max_;        
                return true; 
            }
            else if (tau ==0 || tau < tau_min_) // check this out...
            {
                tau  = tau_min_;
                return true; 
            }
            else{
                return false; 
            }
        }


        bool residual_monotonicity_test(const Vector & g_trial, const Vector & g_old)
        {
            return (norm_l2(g_trial) < norm_l2(g_old)) ? true : false; 
        }


    private:
        VerbosityLevel verbosity_level_;   // verbosity level 

        Scalar tau_max_;            // clamping values of tau to prevent infty 
        Scalar tau_min_;            // clamping values of tau to prevent devision by zero 
        Scalar tau_zero_user_; 

        Scalar alpha_treshold_;     // treshold on scaling
        SizeType max_inner_it_; 

        Matrix I_; 
        bool reset_mass_matrix_; 
        bool is_identity_; 

        Matrix D_; 
        Matrix D_inv_; 
        bool scaling_user_provided_; 
        bool scaling_; 

    };

}
#endif //UTOPIA_ASTRUM_HPP
