#ifndef UTOPIA_AFFINE_SIMILARITY_TRANSFORM_HPP
#define UTOPIA_AFFINE_SIMILARITY_TRANSFORM_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{

    enum AFAlgoVersion  {   AF_VERSION_A  = 1,
                            AF_VERSION_B  = 2,
                            AF_VERSION_C  = 3,
                            AF_VERSION_D  = 4};

    
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class AffineSimilarity : public NonLinearSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef typename NonLinearSolver<Matrix, Vector>::Solver Solver;
        typedef utopia::LSStrategy<Matrix, Vector> LSStrategy; 

    public:
       AffineSimilarity(    const std::shared_ptr <Solver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >(), 
                            const Parameters params                       = Parameters() ):
                            NonLinearSolver<Matrix, Vector>(linear_solver, params), 
                            mass_init_(false), 
                            algo_version_(AF_VERSION_D), 
                            tau_max_(1e9),
                            tau_min_(-1e9), 
                            fix_point_diff_tol_(1e-3),
                            fix_point_max_it_(1000)
                            {
                                //set_parameters(params);

                                verbosity_level_ = params.verbose() ? VERBOSITY_LEVEL_NORMAL : VERBOSITY_LEVEL_QUIET;  
                            }

        bool solve(Function<Matrix, Vector> &fun, Vector &x) override
        {
           using namespace utopia;

           if(mass_init_ == false && mpi_world_rank() == 0)
           {
                std::cerr<<"Affine similarity solver requires mass matrix to be initialized .... \n "; 
                return false; 
           }

            Vector g, s, rhs, x_trial, g_trial;
            Matrix H, A;

            Scalar g_norm, s_norm=9e9, tau;
            SizeType it = 0, it_inner = 0;

            bool converged = false, taken=1;

            fun.gradient(x, g);
            g = -1.0*g;
            g_norm = norm2(g);


            // initialization of  tau 
            tau = 1.0/g_norm; 

            if(verbosity_level_ >= VERBOSITY_LEVEL_NORMAL)
            {
                this->init_solver("Affine similarity", {" it. ", "|| g ||", "|| s_k || ", "tau", "it_inner",  "marker"});
                PrintInfo::print_iter_status(it, {g_norm, 0, tau});
            }

            print_statistics(it, g_norm, tau,  it_inner); 

            it++;

            while(!converged)
            {
                fun.hessian(x, H);
                H *= -1.0; 

                A = H - 1.0/tau * M_; 
                rhs = -1.0 * g; 
                
                //find direction step
                s = local_zeros(local_size(x));
                this->linear_solve(A, rhs, s);

                // build trial point 
                x_trial = x+s; 

                // plain correction, without step-size
                s = 1.0/tau * s; 
                s_norm = norm2(s); 
            
                // gradient of x_trial 
                fun.gradient(x_trial, g_trial);  
                g_trial *= -1.0;     

                if(norm2(g_trial) < norm2(g) && algo_version_ < AF_VERSION_D)
                {
                    x = x_trial; 
                    g = g_trial; 
                    taken = 1; 
                    it_inner = 0; 

                    if(algo_version_ == AF_VERSION_A || algo_version_ == AF_VERSION_C)
                    {
                        tau = estimate_tau(g_trial, g, s, tau, s_norm); 
                        clamp_tau(tau); 
                    }
                }
                else
                {
                    taken=0;
                    bool converged_inner = false; 

                    // here initial value for tau comes from tau_opt 
                    tau = estimate_tau(g_trial, g, s, tau, s_norm); 
                    clamp_tau(tau); 

                    if(verbosity_level_ > VERBOSITY_LEVEL_NORMAL)
                        this->init_solver("Inner it ", {" it. ", "|| tau ||"});
                    
                    if(algo_version_ ==  AF_VERSION_A)
                        converged_inner = true; 

                    Scalar tau_old = 9e9; 
                    it_inner = 0; 

                    while(!converged_inner)
                    {
                        tau_old = tau; 
                        A = H - 1.0/tau * M_; 
                        rhs = -1.0 * g; 
                        
                        //find direction step
                        s = local_zeros(local_size(x));
                        this->linear_solve(A, rhs, s);

                        x_trial = x+s; 

                        // plain correction, without step-size
                        s = 1.0/tau * s; 
                        s_norm = norm2(s); 
                    
                        // gradient of x_trial 
                        fun.gradient(x_trial, g_trial);  
                        g_trial *= -1.0;     

                        tau = estimate_tau(g_trial, g, s, tau, s_norm); 
                        converged_inner =  clamp_tau(tau); 

                        // convergence criterium for fixed point iteration 
                        if(std::abs(tau_old - tau) < fix_point_diff_tol_ || it_inner > fix_point_max_it_ )
                            converged_inner = true; 

                        it_inner++; 

                        if(verbosity_level_ > VERBOSITY_LEVEL_NORMAL)
                            PrintInfo::print_iter_status(it_inner, {tau});
                    }

                    if(norm2(g_trial) < norm2(g))
                    {
                        x = x_trial; 
                        g = g_trial; 

                        // reset value of tau ... 
                        tau = 1.0/g_norm; 
                    }
                    else
                        std::cout<<"--- WARNING: residual monotonicity test failed... \n"; 

                }  // this is outer loop of residual monicity test
            
                g_norm = norm2(g);

                // print iteration status on every iteration
                if(verbosity_level_ >= VERBOSITY_LEVEL_NORMAL)
                    PrintInfo::print_iter_status(it, {g_norm, s_norm, tau, Scalar(it_inner),  Scalar(taken)});

                print_statistics(it, g_norm, tau,  it_inner); 

                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 9e9, s_norm);
                it++;

            } // outer solve loop while(!converged)

            return true;
        }


        void set_mass_matrix(const Matrix & M)
        {
            M_ = M; 
            mass_init_ = true; 
        }


        AFAlgoVersion algo_version() const {  return algo_version_;  }
        void algo_version(const AFAlgoVersion & version) {  algo_version_ = version;  }


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
        Scalar fixed_point_diff_tol() const { return fix_point_diff_tol_; }
        SizeType fixed_point_max_it() const {fix_point_max_it_; }


        void tau_max(const Scalar & tau_max)  {  tau_max_ = tau_max; }
        void tau_min(const Scalar & tau_min)  { tau_min_ = tau_min; }
        void fixed_point_diff_tol(const Scalar & diff_tol)  { fix_point_diff_tol_ = diff_tol; }
        void fixed_point_max_it(const SizeType & max_it_tol)  {fix_point_max_it_ = max_it_tol; }


    private: 
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


        Scalar estimate_tau(const Vector & g_trial, const Vector & g, const Vector & s, const Scalar & tau, const Scalar & s_norm)
        {
            Vector gs_diff = (g_trial - (M_ * s)); 
            Scalar nom = dot(s, ( (1.0/tau * M_ * s) - g)); 
            Scalar help_denom = (2.0 * norm2(gs_diff) * s_norm); 
            return (tau  *  std::abs(nom)/ help_denom); 
        }



        bool clamp_tau(Scalar & tau)
        {
            if(std::isinf(tau) || tau > tau_max_ )
            {
                tau = tau_max_;        
                return true; 
            }
            else if (std::isnan(tau) || tau ==0 || tau < tau_min_)
            {
                tau = tau_max_;        
                return true; 
            }
            else
                return false; 
        }



    private:
        Matrix M_;                  // mass matrix 
        bool mass_init_;            // marker of initialization of mass matrix 

        VerbosityLevel verbosity_level_;   // verbosity level 
        AFAlgoVersion algo_version_; // version 

        Scalar tau_max_;            // clamping values of tau to prevent infty 
        Scalar tau_min_;            // clamping values of tau to prevent devision by zero 
        Scalar fix_point_diff_tol_; // stopping tolerance for fixed point iteration
        SizeType fix_point_max_it_; // maximum iterations for fixed point iteration

    };

}
#endif //UTOPIA_AFFINE_SIMILARITY_TRANSFORM_HPP
