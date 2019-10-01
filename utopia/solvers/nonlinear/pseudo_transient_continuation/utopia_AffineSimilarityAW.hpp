#ifndef UTOPIA_AFFINE_SIMILARITY_AW_HPP
#define UTOPIA_AFFINE_SIMILARITY_AW_HPP

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
    // TODO:: add scaling matrices into evaluations 
    template<class Matrix, class Vector>
    class ODEFormFunction final: public Function<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix);
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::Function<Matrix, Vector>    Fun;

        ODEFormFunction(const std::shared_ptr<Function<Matrix, Vector> > & fun):
        fun_(fun)
        {

        }

        bool value(const Vector &/*point*/, Scalar & value) const override
        {
            value = 0; 
            return false; // should not be necessary
            // TODO:: think about using Normal eq. function
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
    class AffineSimilarityAW : public NewtonBase<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        typedef typename NewtonBase<Matrix, Vector>::Solver Solver;
        typedef utopia::LSStrategy<Vector>                  LSStrategy; 

        using NewtonBase<Matrix, Vector>::print_statistics; 


    public:
       AffineSimilarityAW(  const std::shared_ptr <Solver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >()):
                            NewtonBase<Matrix, Vector>(linear_solver), 
                            tau_max_(1e9),
                            tau_min_(1e-9), 
                            alpha_treshold_(1e-10), 
                            max_inner_it_(5)
                            {
                                verbosity_level_ =  VERBOSITY_LEVEL_NORMAL; 
                            }



        void read(Input &in) override
        {
            NewtonBase<Matrix, Vector>::read(in);
            in.get("tau_max", tau_max_);
            in.get("tau_min", tau_min_);
            in.get("alpha_treshold", alpha_treshold_);
            in.get("max_inner_it", max_inner_it_);
        }

        void print_usage(std::ostream &os) const override
        {
            NewtonBase<Matrix, Vector>::print_usage(os); 
            
        }


        bool solve(Function<Matrix, Vector> &fun_grad, Vector &x) override
        {
           using namespace utopia;

            std::shared_ptr<Function<Matrix, Vector> > fun_grad_ptr_(&fun_grad, [](Function<Matrix, Vector>*){});
            ODEFormFunction<Matrix, Vector> fun(fun_grad_ptr_); 

            Scalar g_norm=0.0, g_norm_old=0.0, s_norm=0.0, mu, L, tau_old; 
            SizeType it_inner = 0; 


            Vector g, s, g_new, r; 
            s = 0 * x; 
            g_new = 0*x; 
            Matrix H; 

            // bool inner_loop = false; 

            fun.hessian(x, H); 

            if(empty(I_)){
                
                if(this->verbose())
                {
                    std::cout<<"mass matrix not set, using Identity matrix ... \n";
                }

                I_  = local_identity(local_size(H).get(0), local_size(H).get(1)); 
            }
            

            if(empty(D_))
            {
                D_ = local_identity(local_size(H).get(0), local_size(H).get(1)); 
                D_inv_ = D_; 
            }


            x = D_inv_ * x; 

            fun.gradient(x, g); 
            g = D_inv_ * g; 
            g_norm = norm2(g); 

            // recompute
            fun.hessian(x, H); 



            Scalar  tau = 1./g_norm; 

            bool converged = false; 
            SizeType it = 0; 

            if(verbosity_level_ >= VERBOSITY_LEVEL_NORMAL)
            {
                this->init_solver("Affine similarity", {" it. ", "|| F ||", "|| Delta x || ", "tau", "mu", "inner it"});
                PrintInfo::print_iter_status(it, {g_norm, 0, tau});
            }

            it++; 

            while(!converged)
            {   
                // to be investigated
                Matrix A = I_ - tau * H; 
                s = 0*x; 
                this->linear_solve(A, g, s);

                r = g - A*s; 
                // std::cout<<"r: "<< norm2(r) << "  \n"; 


                Vector x_trial = x + tau * s; 
                s_norm = norm2(s); 

                fun.gradient(x_trial, g_new); 
                g_new = D_inv_ * g_new; 

                g_norm_old = g_norm; 
                g_norm = norm2(g_new); 


                // tau = estimate_tau(g_new, g, s, tau, s_norm); 

               // Scalar theta = norm_l2(g_new)/norm_l2(g); 


                // if(norm2(g_new) < norm2(g))
                // {
                    mu = estimate_mu(g_new, g, s, tau, s_norm); 
                    // mu = estimate_mu_inexact(g_new, g, s, r, tau, s_norm);

                    // if(mu < 0)
                    // {
                    //     std::cout<<"mu estimated .... \n"; 
                    //     tau = 1./(2. * g_norm_old); 
                    // }
                    // else
                    // {
                        tau = estimate_tau(g_new, g, s, tau, s_norm); 
                        // tau = estimate_tau_inexact(g_new, g, s, r,  tau, s_norm); 
                    // }





                    it_inner = 0; 

                    x = x_trial; 
                    g = g_new; 
                    fun.hessian(x, H); 
                // }
                // else
                // {
                //     mu = estimate_mu(g_new, g, s, tau, s_norm); 
                //     // mu = estimate_mu_inexact(g_new, g, s, r, tau, s_norm);

                //     // if(std::abs(tau-tau_old) < 1e-3){
                //     //     tau = 1./Scalar(norm2(g)); 
                //     // }
                //     // else
                //     // {


                //     // if(mu < 0)
                //     // {
                //     //     std::cout<<"mu estimated2 .... \n"; 
                //     //     tau = 1./(2. * g_norm_old) ; 
                //     // }
                //     // else
                //     // {
                //         tau = estimate_tau(g_new, g, s, tau, s_norm); 
                //         tau_old = tau; 
                //         g_norm = g_norm_old; 
                //         // tau = estimate_tau_inexact(g_new, g, s, r,  tau, s_norm); 
                //     // }

                //     // }

                    
                //     // bool converged_inner = false; 
                //     it_inner++; 
                // }


                // print iteration status on every iteration
                if(verbosity_level_ >= VERBOSITY_LEVEL_NORMAL){
                    PrintInfo::print_iter_status(it, {g_norm, s_norm, tau, mu, Scalar(it_inner)});
                }


                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 9e9, s_norm);
                it++;
            }


            x = D_ * x; 

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

        void tau_max(const Scalar & tau_max)  {  tau_max_ = tau_max; }
        void tau_min(const Scalar & tau_min)  { tau_min_ = tau_min; }

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

    
    public:
        void set_mass_matrix(const Matrix & M)
        {
            I_ = M; 
        }

        void set_scaling_matrix(const Matrix & D)
        {
            D_ = D; 

            Vector d = diag(D_); 
            D_inv_ = diag(1.0/d); 
        }


    private: 
        Scalar estimate_tau(const Vector & g_new, const Vector & g, const Vector & s, const Scalar & tau, const Scalar & /*s_norm*/)
        {   
            Scalar s_norm2 = norm_l2_2(s);

            Scalar nom = dot(g, s) - s_norm2; 
            Scalar denom = dot(g_new, s) - s_norm2; 

            Scalar tau_new = 0.5*tau * std::abs(nom/denom); 
            // bool flg = this->clamp_tau(tau_new); 
            this->clamp_tau(tau_new); 

            if(tau_new==tau_max_ && tau_new == tau)
            {
                tau_new = 1./norm2(g); 
            }
        
            return tau_new; 
        }

        Scalar estimate_tau_inexact(const Vector & g_new, const Vector & g, const Vector & s, const Vector & r,  const Scalar & tau, const Scalar & /*s_norm*/)
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


        Scalar estimate_mu(const Vector & /*g_new*/, const Vector & g, const Vector & s, const Scalar & tau, const Scalar & /*s_norm*/)
        {   
            Scalar s_norm2 = norm_l2_2(s);

            Scalar nom = dot(g, s) - s_norm2; 
            Scalar denom = tau * s_norm2;

            return nom/denom; 
        }


        Scalar estimate_mu_inexact(const Vector & /*g_new*/, const Vector & g, const Vector & s,  const Vector & r, const Scalar & tau, const Scalar & /*s_norm*/)
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

        Scalar alpha_treshold_;     // treshold on scaling
        SizeType max_inner_it_; 

        Matrix I_; 
        Matrix D_; 
        Matrix D_inv_; 

    };

}
#endif //UTOPIA_AFFINE_SIMILARITY_AW_HPP
