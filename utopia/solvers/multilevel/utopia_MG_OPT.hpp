#ifndef UTOPIA_MG_OPT_HPP
#define UTOPIA_MG_OPT_HPP
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NewtonBase.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"
#include "utopia_LS_Strategy.hpp"
#include "utopia_SimpleBacktracking.hpp"

namespace utopia 
{
    /**
     * @brief      The class for Line-search multilevel optimization algorithm. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector>
    class MG_OPT final: public NonlinearMultiLevelBase<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        
        typedef utopia::LSStrategy<Vector>                  LSStrategy; 
        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        
        typedef utopia::NewtonBase<Matrix, Vector>    Solver;
        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;

    

    public:

        bool solve(Vector & x_h) override
        {
            utopia_error("MG_OPT:: solve(x) function is not supported, use solve(fun, x, rhs) instead .... \n"); 
            return false; 
        }

        MG_OPT( const SizeType & n_levels,
                const std::shared_ptr<Solver> &smoother, 
                const std::shared_ptr<Solver> &coarse_solver,
                const std::shared_ptr<LSStrategy> &ls_strategy = std::make_shared<utopia::SimpleBacktracking<Vector> >()): 
                NonlinearMultiLevelBase<Matrix,Vector>(n_levels), 
                _smoother(smoother), 
                _coarse_solver(coarse_solver), 
                _ls_strategy(ls_strategy) 
        {
            
        }
        
        std::string name() override
        {
            return "MG_OPT"; 
        }


        void init_memory(const SizeType & fine_local_size) override 
        {
            std::cout<<"-------- to be done \n"; 
        }


        void read(Input &in) override
        {
            NonlinearMultiLevelBase<Matrix, Vector>::read(in);

            if(_smoother) 
            {
                in.get("smoother", *_smoother);
            }
            if(_coarse_solver) 
            {
                in.get("coarse_solver", *_coarse_solver);
            }     
            if(_ls_strategy) 
            {
                in.get("ls_strategy", *_ls_strategy);
            }                        
        }

        void print_usage(std::ostream &os) const override
        {
            NonlinearMultiLevelBase<Matrix, Vector>::print_usage(os);

            this->print_param_usage(os, "coarse_solver", "NewtonBase", "Input parameters for coarse level QP solvers.", "-"); 
            this->print_param_usage(os, "smoother", "NewtonBase", "Input parameters for fine level QP solver.", "-"); 
            this->print_param_usage(os, "ls_strategy", "LSStrategy", "Input parameters for line-search strategy.", "-"); 
        }  



        bool solve(Fun &fine_fun, Vector & x_h, const Vector & rhs)
        {
            
            bool converged = false;
            SizeType it = 0, n_levels = this->n_levels();
            Scalar r_norm, r0_norm=1, rel_norm=1, energy;

            
            std::string header_message = this->name() + ": " + std::to_string(n_levels) +  " levels";
            this->init_solver(header_message, {" it. ", "|| grad ||", "r_norm" , "Energy"});
            
            this->init_memory(local_size(x_h).get(0)); 
            
            Vector g = local_zeros(local_size(x_h));
            fine_fun.gradient(x_h, g);
            r0_norm = norm2(g);
            r_norm = r0_norm;
            
            fine_fun.value(x_h, energy);
            
            if(this->verbose())
                PrintInfo::print_iter_status(it, {r_norm, rel_norm, energy});
            
            it++;
            
            while(!converged)
            {

                this->multiplicative_cycle(fine_fun, x_h, rhs, n_levels);
                
#ifdef CHECK_NUM_PRECISION_mode
                if(has_nan_or_inf(x_h) == 1)
                {
                    x_h = local_zeros(local_size(x_h));
                    return true;
                }
#endif
                
                fine_fun.gradient(x_h, g);
                fine_fun.value(x_h, energy);
                
                r_norm = norm2(g);
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
            if(has_nan_or_inf(x_h) == 1)
                exit(0);
#endif
            
            
            return true;
        }



    private: 

        bool multiplicative_cycle(Fun &fine_fun, Vector & u_l, const Vector &f, const SizeType & l)
        {
           Vector g_fine, g_restricted,  g_coarse, u_2l, s_coarse, s_fine, u_init; 
           Scalar alpha; 

           this->make_iterate_feasible(fine_fun, u_l); 

            // PRE-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->pre_smoothing_steps()); 
            fine_fun.gradient(u_l, g_fine);             

            g_fine -= f; 

            this->transfer(l-2).restrict(g_fine, g_restricted); 
            this->transfer(l-2).project_down(u_l, u_2l); 

            this->make_iterate_feasible(this->function(l-2), u_2l); 
            this->zero_correction_related_to_equality_constrain(this->function(l-2), g_restricted); 

            this->function(l-2).gradient(u_2l, g_coarse); 

            u_init = u_2l; 
            g_coarse -= g_restricted;  // tau correction - g_diff in rmtr 
                                 

            if(l == 2)
            {
                coarse_solve(this->function(0), u_2l, g_coarse); 
            }
            else
            {
                // recursive call into FAS - needs to be checked 
                for(SizeType k = 0; k < this->mg_type(); k++)
                {   
                    SizeType l_new = l - 1; 
                    this->multiplicative_cycle(this->function(l-2), u_2l, g_coarse, l_new); 
                }
            }

            s_coarse = u_2l - u_init; 
            this->transfer(l-2).interpolate(s_coarse, s_fine);
            this->zero_correction_related_to_equality_constrain(fine_fun, s_fine); 
            

            // TODO:: check which gradient: with/without tau correction  ...
            // TODO::  put  energy checks ... 


            _ls_strategy->get_alpha(fine_fun, g_fine, u_l, s_fine, alpha);
            u_l += alpha * s_fine; 

            // POST-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->post_smoothing_steps()); 

            return true; 
        }


        /**
         * @brief      The function invokes smoothing. 
         *
         * @param[in]  A     The stifness matrix. 
         * @param[in]  rhs   The right hand side.
         * @param      x     The current iterate. 
         * @param[in]  nu    The number of smoothing steps. 
         *
         * @return   
         */
        bool smoothing(Fun &fun,  Vector &x, const Vector &rhs, const SizeType & nu = 1)
        {
            _smoother->max_it(1); 
            _smoother->verbose(true); 
            _smoother->solve(fun, x, rhs); 

            return true; 
        }

        /**
         * @brief      The function invokes coarse solver. 
         *
         * @param      fun   Function to be minimized 
         * @param[in]  rhs   The right hand side.
         * @param      x     The current iterate. 
         * @param[in]  nu    The number of smoothing steps. 
         *
         * @return   
         */
        bool coarse_solve(Fun &fun, Vector &x, const Vector & rhs)
        {
            _coarse_solver->verbose(true); 
            _coarse_solver->solve(fun, x, rhs); 
            return true; 
        }


    public:
        /**
         * @brief      Function changes direct solver needed for coarse grid solve. 
         *
         * @param[in]  linear_solver  The linear solver.
         *
         * @return     
         */
        bool set_coarse_solver(const std::shared_ptr<Solver> &nonlinear_solver = std::shared_ptr<Solver>())
        {
            _coarse_solver = nonlinear_solver; 
            // _coarse_solver->set_parameters(_parameters); 
            return true; 
        }

        /**
         * @brief      Function changes direct solver needed for coarse grid solve. 
         *
         * @param[in]  linear_solver  The linear solver.
         *
         * @return     
         */
        bool set_smoother(const std::shared_ptr<Solver> &nonlinear_solver = std::shared_ptr<Solver>())
        {
            _smoother = nonlinear_solver; 
            // _smoother->set_parameters(_parameters); 
            return true; 
        }


        /**
         * @brief      Sets line search strategy. 
         */
        bool set_line_search(const std::shared_ptr<LSStrategy> &ls_strategy = std::shared_ptr<LSStrategy>())
        {
            _ls_strategy = ls_strategy; 
            return true; 
        }


    protected:   
        std::shared_ptr<Solver>             _smoother;          /*  optimization method used on fine levels */
        std::shared_ptr<Solver>             _coarse_solver;     /*  optimization method for coarse level */
        std::shared_ptr<LSStrategy>         _ls_strategy;       /*  LS used to determine step size inside of MG */



    };

}

#endif //UTOPIA_MG_OPT_HPP

