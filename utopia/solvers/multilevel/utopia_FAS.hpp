#ifndef UTOPIA_FAS_HPP
#define UTOPIA_FAS_HPP
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"

#include "utopia_LevelMemory.hpp"


namespace utopia 
{
    /**
     * @brief      The class for Full approximation scheme.
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector>
    class FAS : public NonlinearMultiLevelBase<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::NonLinearSolver<Matrix, Vector>     Solver;
        typedef utopia::NonLinearSmoother<Matrix, Vector>   Smoother;
        typedef utopia::Transfer<Matrix, Vector>   Transfer;
        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;

        public:

            FAS( const std::shared_ptr<Smoother> &smoother, const std::shared_ptr<Solver> &coarse_solver, const Parameters params = Parameters()): 
                    NonlinearMultiLevelBase<Matrix,Vector>(params),
                    smoother_(smoother), 
                    coarse_solver_(coarse_solver) 
            {
                set_parameters(params); 
            }

            virtual ~FAS(){} 
        

            void set_parameters(const Parameters params) override
            {
                NonlinearMultiLevelBase<Matrix, Vector>::set_parameters(params); 
                smoother_->set_parameters(params); 
                coarse_solver_->set_parameters(params); 
                
                parameters_ = params; 
            }


            virtual std::string name_id() override
            {
                return "FAS"; 
            }

    
        virtual bool solve(Fun &fine_fun, Vector & x_h, const Vector & rhs) override
        {

            std::cout<<"yes, this solve ------- \n";
            
            bool converged = false;
            SizeType it = 0, n_levels = this->n_levels();
            Scalar r_norm, r0_norm=1, rel_norm=1, energy;

            
            std::string header_message = this->name_id() + ": " + std::to_string(n_levels) +  " levels";
            this->init_solver(header_message, {" it. ", "|| grad ||", "r_norm" , "Energy"});
            
            this->status_.clear();

            SizeType loc_size = local_size(x_h).get(0);
            this->init_memory(loc_size); 


            memory_.x[n_levels-1] = x_h;
            memory_.g[n_levels-1] = local_zeros(local_size(memory_.x[n_levels-1]));


            fine_fun.gradient(memory_.x[n_levels-1], memory_.g[n_levels-1]);
            r0_norm = norm2(memory_.g[n_levels-1]);
            r_norm = r0_norm;
            
            fine_fun.value(memory_.x[n_levels-1], energy);
            
            if(this->verbose())
                PrintInfo::print_iter_status(it, {r_norm, rel_norm, energy});
            
            it++;
            
            while(!converged)
            {
                this->multiplicative_cycle(fine_fun, memory_.x[n_levels-1], rhs, n_levels);
                
#ifdef CHECK_NUM_PRECISION_mode
                if(has_nan_or_inf(memory_.x[n_levels-1]) == 1)
                {
                    memory_.x[n_levels-1] = local_zeros(local_size(memory_.x[n_levels-1]));
                    return true;
                }
#endif
                
                fine_fun.gradient(memory_.x[n_levels-1], memory_.g[n_levels-1]);
                fine_fun.value(memory_.x[n_levels-1], energy);
                
                r_norm = norm2(memory_.g[n_levels-1]);
                rel_norm = r_norm/r0_norm;
                
                // print iteration status on every iteration
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {r_norm, rel_norm, energy});
                
                // check convergence and print interation info
                converged = this->check_convergence(it, r_norm, rel_norm, 1);
                it++;
            }
            
            this->print_statistics(it);
            
#ifdef CHECK_NUM_PRECISION_mode
            if(has_nan_or_inf(memory_.x[n_levels-1]) == 1)
                exit(0);
#endif
            x_h = memory_.x[n_levels-1]; 
            return true;
        }



    protected: 

        virtual void init_memory(const SizeType & fine_local_size) override 
        {
            memory_.init(this->n_levels()); 
            memory_.g_diff[this->n_levels()-1] = local_zeros(fine_local_size); 
        }


        bool multiplicative_cycle(Fun &fine_fun, Vector & u_l, const Vector &f, const SizeType & l) override
        {
            for(auto l = this->n_levels()-1; l > 0; l--)
            {
                // pre-smoothing 
                smoothing(this->function(l), memory_.x[l], memory_.g_diff[l], this->pre_smoothing_steps()); 

                this->transfer(l-1).project_down(memory_.x[l], memory_.x[l-1]); 
                memory_.x_0[l-1] = memory_.x[l-1]; 


                // multilevel gradient ... 
                this->function(l).gradient(memory_.x[l], memory_.g[l]); 
                memory_.g[l] -= memory_.g_diff[l]; 


                this->transfer(l-1).restrict(memory_.g[l], memory_.g_diff[l-1]);

                this->function(l-1).gradient(memory_.x[l-1], memory_.g[l-1]); 
                memory_.g_diff[l-1] = memory_.g[l-1] - memory_.g_diff[l-1]; 

                this->zero_correction_related_to_equality_constrain(this->function(l-1), memory_.g_diff[l-1]); 
            }

            coarse_solve(this->function(0), memory_.x[0],memory_.g_diff[0]); 

            for(auto l = 0; l < this->n_levels()-1; l++)
            {
                memory_.c[l] = memory_.x[l] - memory_.x_0[l]; 
                this->transfer(l).interpolate(memory_.c[l], memory_.c[l+1]);

                this->zero_correction_related_to_equality_constrain(this->function(l+1), memory_.c[l+1]); 

                memory_.x[l+1] += memory_.c[l+1]; 
                smoothing(this->function(l+1), memory_.x[l+1], memory_.g_diff[l+1], this->pre_smoothing_steps()); 
            }

            return true; 
        }

        bool smoothing(Function<Matrix, Vector> &fun,  Vector &x, const Vector &f, const SizeType & nu = 1)
        {
            smoother_->sweeps(nu); 
            smoother_->nonlinear_smooth(fun, x, f); 
            return true; 
        }

        bool coarse_solve(Fun &fun, Vector &x, const Vector & rhs) override
        {   
            coarse_solver_->max_it(1); 
            coarse_solver_->solve(fun, x, rhs); 
            return true; 
        }

    public:
        bool change_coarse_solver(const std::shared_ptr<Solver> &nonlinear_solver = std::shared_ptr<Solver>())
        {
            coarse_solver_ = nonlinear_solver; 
            return true; 
        }



    protected:   
        std::shared_ptr<Smoother>           smoother_;
        std::shared_ptr<Solver>             coarse_solver_;  
        LevelMemory <Matrix, Vector>         memory_;

    private:
        Parameters                          parameters_; 


    };

}

#endif //UTOPIA_FAS_HPP

