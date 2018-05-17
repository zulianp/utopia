#ifndef UTOPIA_FAS_HPP
#define UTOPIA_FAS_HPP
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"


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
                    NonlinearMultiLevelBase<Matrix,Vector>(params), _smoother(smoother), _coarse_solver(coarse_solver) 
            {
                set_parameters(params); 
            }

            virtual ~FAS(){} 
        

            void set_parameters(const Parameters params) override
            {
                NonlinearMultiLevelBase<Matrix, Vector>::set_parameters(params); 
                _smoother->set_parameters(params); 
                _coarse_solver->set_parameters(params); 
                
                _parameters = params; 
            }


            virtual std::string name_id() override
            {
                return "FAS"; 
            }

    

    private: 



        bool multiplicative_cycle(Fun &fine_fun, Vector & u_l, const Vector &f, const SizeType & l) override
        {
            this->memory().x[this->n_levels()-1] = u_l; 

            // sto be investigated with the energy  ... 
            // this->make_iterate_feasible(this->function(this->n_levels()-1), memory.x[this->n_levels()-1]); 

            for(auto l = this->n_levels()-1; l > 0; l--)
            {
                smoothing(this->function(l), this->memory().x[l], this->memory().g_diff[l], this->pre_smoothing_steps()); 

                this->transfer(l-1).project_down(this->memory().x[l], this->memory().x[l-1]); 
                this->memory().x_0[l-1] = this->memory().x[l-1]; 

                // TODO:: make this function nicer... 
                // generic enough...
                this->get_multilevel_gradient(l); 



                this->transfer(l-1).restrict(this->memory().g[l], this->memory().g_diff[l-1]);

                this->function(l-1).gradient(this->memory().x[l-1], this->memory().g[l-1]); 
                this->memory().g_diff[l-1] = this->memory().g[l-1] - this->memory().g_diff[l-1]; 

                this->zero_correction_related_to_equality_constrain(this->function(l-1), this->memory().g_diff[l-1]); 
            }

            coarse_solve(this->function(0), this->memory().x[0], this->memory().g_diff[0]); 

            for(auto l = 0; l < this->n_levels()-1; l++)
            {
                this->memory().c[l] = this->memory().x[l] - this->memory().x_0[l]; 
                this->transfer(l).interpolate(this->memory().c[l], this->memory().c[l+1]);

                this->zero_correction_related_to_equality_constrain(this->function(l+1), this->memory().c[l+1]); 

                this->memory().x[l+1] += this->memory().c[l+1]; 
                smoothing(this->function(l+1), this->memory().x[l+1], this->memory().g_diff[l+1], this->pre_smoothing_steps()); 
            }


            // to be fixed...
            u_l = this->memory().x[this->n_levels()-1]; 

            return true; 
        }



        // TODO:: find more suitable place for  this function 
        // ideally, such that it also works for 
        virtual bool get_multilevel_gradient(const SizeType & level)
        {
            if(level < this->n_levels())
            {
                this->function(level).gradient(this->memory().x[level], this->memory().g[level]); 
                this->memory().g[level] -= this->memory().g_diff[level]; 

                return true; 
            }
            else
            {
                return this->function(level).gradient(this->memory().x[level], this->memory().g[level]); 
            }
        }



        bool smoothing(Function<Matrix, Vector> &fun,  Vector &x, const Vector &f, const SizeType & nu = 1)
        {
            _smoother->sweeps(nu); 
            _smoother->nonlinear_smooth(fun, x, f); 
            return true; 
        }

        bool coarse_solve(Fun &fun, Vector &x, const Vector & rhs) override
        {   
            _coarse_solver->max_it(1); 
            _coarse_solver->solve(fun, x, rhs); 
            return true; 
        }

    public:
        bool change_coarse_solver(const std::shared_ptr<Solver> &nonlinear_solver = std::shared_ptr<Solver>())
        {
            _coarse_solver = nonlinear_solver; 
            return true; 
        }



    protected:   
        std::shared_ptr<Smoother>           _smoother;
        std::shared_ptr<Solver>             _coarse_solver;  

    private:
        Parameters                          _parameters; 


    };

}

#endif //UTOPIA_FAS_HPP

