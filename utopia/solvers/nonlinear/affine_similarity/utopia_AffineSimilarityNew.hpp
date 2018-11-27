#ifndef UTOPIA_AFFINE_SIMILARITY_NEW_HPP
#define UTOPIA_AFFINE_SIMILARITY_NEW_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_LS_Strategy.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{
    
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class AffineSimilarityNew : public NewtonBase<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;
        typedef typename NewtonBase<Matrix, Vector>::Solver Solver;
        typedef utopia::LSStrategy<Vector>                  LSStrategy; 

    public:
       AffineSimilarityNew( const std::shared_ptr <Solver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >(), 
                            const Parameters params                       = Parameters() ):
                            NewtonBase<Matrix, Vector>(linear_solver, params), 
                            mass_init_(false), 
                            tau_max_(1e14),
                            tau_min_(1e-14)
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

            bool converged = false; 

            gradient(fun, x, g); 
            g_norm = norm2(g);

            SizeType solves_counter = 0; 

            // initialization of  tau 
            tau = 1.0/g_norm;  

            if(verbosity_level_ >= VERBOSITY_LEVEL_NORMAL)
            {
                this->init_solver("Affine similarity", {" it. ", "|| F ||", "|| Delta x || ", "tau", "it_inner"});
                PrintInfo::print_iter_status(it, {g_norm, 0, tau});
            }

            print_statistics(it, g_norm, tau,  it_inner); 
            it++;

            while(!converged)
            {
                hessian(fun, x, H); 
                gradient(fun, x, g); 

                A = 1./tau* M_ - H; 
                rhs = g; 

                //find direction step
                s = local_zeros(local_size(x));
                this->linear_solve(A, rhs, s);
                solves_counter++; 

                // build trial point 
                x_trial = x + s; 

                // correction without time step 
                s = 1./tau * s; 
                s_norm = norm2(s); 

                // gradient of trial point 
                gradient(fun, x_trial, g_trial);                 

                if(norm2(g_trial) < norm2(g))
                {   
                    // update tau since, you  have already all ingredients 
                    tau = estimate_tau(g_trial, g, s, tau, s_norm); 
                    clamp_tau(tau); 

                    x = x_trial; 
                    it_inner = 1; 
                }
                else
                { 
                    if(verbosity_level_ > VERBOSITY_LEVEL_NORMAL)
                    {
                        this->init_solver("Fixed point it ", {" it. ", "|| tau ||"});
                        PrintInfo::print_iter_status(0, {tau});
                    }
                    
                    bool converged_inner = false;
                    it_inner = 1; 

                    if(verbosity_level_ > VERBOSITY_LEVEL_NORMAL)
                        PrintInfo::print_iter_status(it_inner, {tau});

                    while(!converged_inner && it_inner < 10)
                    {
                        // here initial value for tau comes from tau_opt 
                        tau = estimate_tau(g_trial, g, s, tau, s_norm); 
                        converged_inner = clamp_tau(tau); 

                        A = 1./tau* M_ - H; 

                        //find direction step
                        s = local_zeros(local_size(x));
                        this->linear_solve(A, rhs, s);
                        solves_counter++; 
                        it_inner++;

                        // build trial point 
                        x_trial = x + s; 

                        // correction without time step 
                        s = 1./tau * s; 
                        s_norm = norm2(s); 

                        // gradient of trial point 
                        gradient(fun, x_trial, g_trial);                            

                        // we iterate until residual monotonicity test is satisfied... 
                        if(norm2(g_trial) < norm2(g))
                        {   
                            x = x_trial; 
                            converged_inner = true; 
                        }

                        if(!converged_inner && verbosity_level_ > VERBOSITY_LEVEL_NORMAL)
                            PrintInfo::print_iter_status(it_inner, {tau});
                    }

                    if (mpi_world_rank() == 0 && verbosity_level_ > VERBOSITY_LEVEL_NORMAL)
                        std::cout<<"------------------------------ end of fixed point iteration ------------------------ \n"; 

                }  // this is outer loop of residual monicity test
                
                g_norm = norm2(g_trial);

                // print iteration status on every iteration
                if(verbosity_level_ >= VERBOSITY_LEVEL_NORMAL)
                    PrintInfo::print_iter_status(it, {g_norm, s_norm, tau, Scalar(it_inner-1)});

                print_statistics(it, g_norm, tau,  it_inner); 

                // check convergence and print interation info
                converged = this->check_convergence(it, g_norm, 9e9, s_norm);
                it++;

            } // outer solve loop while(!converged)

            if (mpi_world_rank() == 0)
                std::cout<<"solves_counter: "<< solves_counter << "  \n"; 

            return true;
        }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        void set_mass_matrix(const Matrix & M)
        {
            M_      = M; 

            Vector d = diag(M_); 
            M_inv_  = diag( 1.0/ d); 
            c_      = max(d); 

            mass_init_ = true; 
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


        virtual void print_statistics(const SizeType & it_global) override
        {
            NonLinearSolver<Matrix, Vector>::print_statistics(it_global); 
        }
    
    private: 
        // Scalar estimate_tau(const Vector & g_trial, const Vector & g, const Vector & s, const Scalar & tau, const Scalar & s_norm)
        // {   
        //     Vector gs_diff = (g_trial - (M_ * s)); 
        //     Scalar nom = dot(s, ( (1.0/tau * M_ * s) - g)); 
        //     Scalar help_denom = (2.0 * norm2(gs_diff) * s_norm); 
        //     return (tau  *  std::abs(nom)/ help_denom); 
        // }


        Scalar estimate_tau(const Vector & g_trial, const Vector & g, const Vector & s, const Scalar & tau, const Scalar & s_norm)
        {   
            Scalar nom = std::abs(dot(s, ( (M_ * s) - g))) * tau; 
            
            Vector gs_diff = ((M_inv_ * g_trial) - s); 
            Scalar help_denom = (2.0 * norm2(gs_diff) * s_norm) * c_; 
            
            return nom/ help_denom; 
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
                tau = tau_min_;        
                return true; 
            }
            else
                return false; 
        }


        void gradient(Function<Matrix, Vector> &fun,  const Vector & x, Vector & g)
        {
            fun.gradient(x, g); 
            g = -1.0 * g; // because we work with negative definitness 
        }

        void hessian(Function<Matrix, Vector> &fun,  const Vector & x, Matrix & H)
        {
            fun.hessian(x, H); 
            H = -1.0 * H; // because we work with negative definitness
        }


    private:
        Matrix M_;                  // mass matrix 
        Matrix M_inv_;                  // inverse of mass matrix 
        Scalar c_; 


        bool mass_init_;            // marker of initialization of mass matrix 

        VerbosityLevel verbosity_level_;   // verbosity level 

        Scalar tau_max_;            // clamping values of tau to prevent infty 
        Scalar tau_min_;            // clamping values of tau to prevent devision by zero 

    };

}
#endif //UTOPIA_AFFINE_SIMILARITY_NEW_HPP
