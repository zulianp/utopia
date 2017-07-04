/*
* @Author: alenakopanicakova
* @Date:   2017-05-03
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-04
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
        
        virtual std::string name_id()
        {
            return "MG_OPT"; 
        }

        void set_parameters(const Parameters params)  // override
        {
            NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::set_parameters(params); 
            _smoother->set_parameters(params); 
            _coarse_solver->set_parameters(params); 
            _parameters = params; 
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


        bool multiplicative_cycle(FunctionType &fine_fun, Vector & u_l, const Vector &f, const SizeType & l) override
        {
           Vector g_fine, g_coarse, u_2l, e, u_init; 
           Scalar alpha; 

           this->make_iterate_feasible(fine_fun, u_l); 

            // PRE-SMOOTHING 
            smoothing(fine_fun, u_l, f, this->pre_smoothing_steps()); 
            fine_fun.gradient(u_l, g_fine);             

            g_fine -= f; 

            transfers(l-2).restrict(g_fine, g_fine); 
            transfers(l-2).project_down(u_l, u_2l); 

            this->make_iterate_feasible(levels(l-2), u_2l); 
            this->zero_correction_related_to_equality_constrain(levels(l-2), g_fine); 

            levels(l-2).gradient(u_2l, g_coarse); 

            u_init = u_2l; 
            g_coarse -= g_fine;  // tau correction 
                                 
            if(l == 2)
            {
                this->make_iterate_feasible(levels(0), u_2l); 
                coarse_solve(levels(0), u_2l, g_coarse); 
            }
            else
            {
                // recursive call into FAS - needs to be checked 
                for(SizeType k = 0; k < this->mg_type(); k++)
                {   
                    SizeType l_new = l - 1; 
                    this->multiplicative_cycle(levels(l-2), u_2l, g_coarse, l_new); 
                }
            }

            e = u_2l - u_init; 
            transfers(l-2).interpolate(e, e);
            this->zero_correction_related_to_equality_constrain(fine_fun, e); 
            
            _ls_strategy->get_alpha(fine_fun, g_fine, u_l, e, alpha);
            // std::cout<<"l:   "<< l << "   alpha:   "<< alpha << " <g_k, e_h>:   "<< dot(L_l, e_h) << "  \n"; 

            u_l += alpha * e; 

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
        bool coarse_solve(FunctionType &fun, Vector &x, const Vector & rhs) override
        {
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
        std::shared_ptr<Solver>             _coarse_solver;     /*  optimization method for coarse level */
        std::shared_ptr<LSStrategy>         _ls_strategy;       /*  LS used to determine step size inside of MG */


    private:
        Parameters                          _parameters; 


    };

}

#endif //UTOPIA_MG_OPT_HPP

