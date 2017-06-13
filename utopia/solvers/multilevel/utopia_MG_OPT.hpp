/*
* @Author: alenakopanicakova
* @Date:   2017-05-03
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-06-13
*/

#ifndef UTOPIA_MG_OPT_HPP
#define UTOPIA_MG_OPT_HPP
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"
#include "utopia_LS_Strategy.hpp"

namespace utopia 
{
    /**
     * @brief      The class for Line-search multilevel optimization algorithm. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector, class FunctionType>
    class MG_OPT : public NonlinearMultiLevelBase<Matrix, Vector, FunctionType>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::NonLinearSolver<Matrix, Vector>     Solver;
        typedef utopia::LSStrategy<Matrix, Vector>          LSStrategy; 
        typedef utopia::Transfer<Matrix, Vector>            Transfer;

    

    public:

        MG_OPT( const std::shared_ptr<Solver> &smoother = std::shared_ptr<Solver>(), 
                const std::shared_ptr<Solver> &coarse_solver = std::shared_ptr<Solver>(),
                const std::shared_ptr<LSStrategy> &ls_strategy = std::shared_ptr<LSStrategy>(),
                const Parameters params = Parameters()): 
                NonlinearMultiLevelBase<Matrix,Vector, FunctionType>(params), 
                _smoother(smoother), 
                _coarse_solver(coarse_solver), 
                _ls_strategy(ls_strategy) 
        {
            set_parameters(params); 
        }

        virtual ~MG_OPT(){} 
        

        void set_parameters(const Parameters params)  // override
        {
            NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::set_parameters(params); 
            _smoother->set_parameters(params); 
            _coarse_solver->set_parameters(params); 
            
            _parameters = params; 
        }

        virtual bool solve(FunctionType & fine_fun, Vector &x_h)
        {
            Vector rhs = local_zeros(local_size(x_h)); 
            return solve(fine_fun,  x_h, rhs); 
        }

        /**
         * @brief      The solve function for multigrid method. 
         *
         * @param[in]  rhs   The right hand side.
         * @param      x_0   The initial guess. 
         *
         */
        virtual bool solve(FunctionType &fine_fun, Vector & x_h, const Vector & rhs) 
        {
            this->init_solver("MG_OPT", {" it. ", "|| r_N ||", "r_norm" , "E"}); 

            Vector F_h  = local_zeros(local_size(x_h)); 

            bool converged = false; 
            SizeType it = 0, l = this->num_levels(); 
            Scalar r_norm, r0_norm, rel_norm;

            std::cout<<"MG_OPT: number of levels: "<< l << "  \n"; 

            // just to check what is problem 
            Matrix hessian; 
            fine_fun.hessian(x_h, hessian); 

            fine_fun.gradient(x_h, F_h); 
            r0_norm = norm2(F_h); 

            if(this->verbose())
                PrintInfo::print_iter_status(it, {r0_norm}); 

            it++; 

            while(!converged)
            {            
                if(this->cycle_type() =="multiplicative")
                    multiplicative_cycle(fine_fun, x_h, rhs, l); 
                // else if(this->cycle_type() =="full")
                //     full_cycle(rhs, l, x_0); 
                else
                    std::cout<<"ERROR::MG_OPT<< unknown MG type... \n"; 


                #ifdef CHECK_NUM_PRECISION_mode
                    if(has_nan_or_inf(x_h) == 1)
                    {
                        x_h = local_zeros(local_size(x_h));
                        return true; 
                    }
                #endif    

                fine_fun.gradient(x_h, F_h); 
                
                Scalar energy; 
                fine_fun.value(x_h, energy); 
                
                r_norm = norm2(F_h);
                rel_norm = r_norm/r0_norm; 

                // print iteration status on every iteration 
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {r_norm, rel_norm, energy}); 

                // check convergence and print interation info
                converged = this->check_convergence(it, r_norm, rel_norm, 1); 
                it++; 
            
            }

            return true; 
        }



    private: 

        inline FunctionType &levels(const SizeType &l)
        {
            return this->_nonlinear_levels[l]; 
        }

        inline Transfer &transfers(const SizeType & l)
        {
            return this->_transfers[l]; 
        }

        bool full_cycle(FunctionType &fine_fun, Vector & u_l, const Vector &f, const SizeType & l)
        {
            std::cout<<"-------- not yet....... \n"; 
            return true; 
        }



        bool multiplicative_cycle(FunctionType &fine_fun, Vector & u_l, const Vector &f, const SizeType & l)
        {
            Vector L_l, L_2l, r_h,  r_2h, u_2l, e_2h, e_h, u_init; 
            Scalar alpha; 

            this->make_iterate_feasible(fine_fun, u_l); 

            // PRE-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->pre_smoothing_steps()); 
            fine_fun.gradient(u_l, L_l);             

            r_h = L_l - f; 

            transfers(l-2).restrict(r_h, r_2h); 
            transfers(l-2).project_down(u_l, u_2l); 

            this->make_iterate_feasible(levels(l-2), u_2l); 
            this->zero_boundary_correction(levels(l-2), r_2h); 

            levels(l-2).gradient(u_2l, L_2l); 

            u_init = u_2l; 
            L_2l = L_2l - r_2h;  // tau correction 
                                 
            if(l == 2)
            {
                this->make_iterate_feasible(levels(0), u_2l); 
                coarse_solve(levels(0), u_2l, L_2l); 
            }
            else
            {
                for(SizeType k = 0; k < this->mg_type(); k++)
                {   
                    SizeType l_new = l - 1; 
                    multiplicative_cycle(levels(l-2), u_2l, L_2l, l_new); 
                }
            }

            e_2h = u_2l - u_init; 
            transfers(l-2).interpolate(e_2h, e_h);

            this->zero_boundary_correction(fine_fun, e_h); 

            fine_fun.gradient(u_l, L_l); 

            _ls_strategy->get_alpha(fine_fun, L_l, u_l, e_h, alpha);
            // std::cout<<"l:   "<< l << "   alpha:   "<< alpha << " <g_k, e_h>:   "<< dot(L_l, e_h) << "  \n"; 


            u_l += alpha * e_h; 

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
        bool smoothing(Function<Matrix, Vector> &fun,  Vector &x, const Vector &rhs, const SizeType & nu = 1)
        {
            _smoother->max_it(nu); 
            // _smoother->verbose(true); 
            _smoother->solve(fun, x, rhs); 
            return true; 
        }


        bool coarse_solve(FunctionType &fun, Vector &x, const Vector & rhs)
        {
            // _coarse_solver->verbose(true); 
            // _coarse_solver->max_it(10); 
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
            _coarse_solver->set_parameters(_parameters); 
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
            _smoother->set_parameters(_parameters); 
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
        std::shared_ptr<Solver>             _coarse_solver;     /*  optimization method on coarse level */
        std::shared_ptr<LSStrategy>         _ls_strategy;       /*  LS used to determine step size inside of MG */


    private:
        Parameters                          _parameters; 


    };

}

#endif //UTOPIA_MG_OPT_HPP

