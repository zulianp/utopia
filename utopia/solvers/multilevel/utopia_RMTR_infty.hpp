/*
* @Author: alenakopanicakova
* @Date:   2017-04-19
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-06-15
*/

#ifndef UTOPIA_RMTR_INFTY_HPP
#define UTOPIA_RMTR_INFTY_HPP
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"

#include "utopia_TRBoxSubproblem.hpp"
#include "utopia_Linear.hpp"
#include "utopia_Level.hpp"
#include "utopia_LS_Strategy.hpp"

#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"

#include "utopia_MultiLevelEvaluations.hpp"

namespace utopia 
{

    /**
     * @brief      The class for Nonlinear Multigrid solver. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector, class FunctionType, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER >
    class RMTR_infty : public RMTR<Matrix, Vector, FunctionType, CONSISTENCY_LEVEL>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;
        typedef utopia::NonLinearSolver<Matrix, Vector>     Solver;
        typedef utopia::NonLinearSmoother<Matrix, Vector>   Smoother;
        typedef utopia::TRBoxSubproblem<Matrix, Vector>        TRBoxSubproblem; 
        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef utopia::Level<Matrix, Vector>               Level;

    
    public:

       /**
        * @brief      Multigrid class
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level. 
        */
        RMTR_infty(    
                const std::shared_ptr<TRBoxSubproblem> &tr_subproblem_coarse = std::shared_ptr<TRBoxSubproblem>(),
                const std::shared_ptr<TRBoxSubproblem> &tr_subproblem_smoother = std::shared_ptr<TRBoxSubproblem>(),
                // const std::shared_ptr<LSStrategy> &ls_strategy = std::shared_ptr<LSStrategy>(),
                const Parameters params = Parameters()): 
                RMTR<Matrix,Vector, FunctionType, CONSISTENCY_LEVEL>(tr_subproblem_coarse, tr_subproblem_smoother, params) 
        {
            set_parameters(params); 
        }

        virtual ~RMTR_infty(){} 
        

        void set_parameters(const Parameters params)  // override
        {
            RMTR<Matrix, Vector, FunctionType, CONSISTENCY_LEVEL>::set_parameters(params);    
        }


        std::string name_id() { return "RMTR_infty";  }
        

    private: 

        inline FunctionType &levels(const SizeType &l)
        {
            return this->_nonlinear_levels[l]; 
        }

        inline Transfer &transfers(const SizeType & l)
        {
            return this->_transfers[l]; 
        }

 protected:
        // virtual void delta_update(const Scalar & rho, const SizeType & level, const Vector & s_global, bool & converged) override
        // {
        //     return; 
        // }



        // virtual bool check_local_convergence(const SizeType & it_success, const Scalar & g_norm, const SizeType & level, const Scalar & delta) override
        // {   
        //     // coarse one 
        //     if(level == 1 && (it_success >= _max_coarse_it))
        //         return true; 
        //     // every other level 
        //     else if (level > 1 && it_success >= _max_smoothing_it)
        //         return true; 

        //     if(delta < _delta_min)
        //         return true; 

        //     return criticality_measure_termination(g_norm); 
        
        //     // return true; 
        // }

        // /**
        //  * @brief      THis check guarantee that iterates at a lower level remain in the TR radius defined at the calling level
        //  *
        //  * @param[in]  corr_norm  The norm of sum of all corrections on given level 
        //  * @param[in]  level      The level
        //  *
        //  */
        // virtual bool delta_termination(const Scalar & corr_norm, const SizeType & level) override
        // {   
        //     return (corr_norm > (1 - _eps_delta_termination) * get_delta(level-1)) ? true : false; 
        //     // return true; 
        // }



        // // needs to be overriden in case of \infty norm 
        // virtual bool criticality_measure_termination(const Scalar & g_norm) override
        // {
        //     // this should be fancier based on Graffon paper, but it is quite boring to work with so many mesh informations
        //     // Scalar _eps_grad_termination = std::min(0.001, eps_grad_termination[level]/ psi); 
            
        //     return (g_norm < _eps_grad_termination) ? true : false;    
            
        //     // return true; 
        // }



        // virtual bool grad_smoothess_termination(const Vector & g_restricted, const Vector & g_coarse) override
        // {
        //     Scalar Rg_norm = norm2(g_restricted); 
        //     Scalar g_norm = norm2(g_coarse);
        //     // return (Rg_norm >= _grad_smoothess_termination g_norm) ? true : false;   
            
        //     if(Rg_norm >= _grad_smoothess_termination * g_norm)
        //     {
        //         std::cout<<"grad is not smooth enough............... => NO RECURSION ANYMORE ....  \n"; 
        //         return true; 
        //     }
        //     else
        //         return false; 


        //     // return true; 

        // }



        // virtual bool solve_qp_subproblem(const Matrix & H, const Vector & g, Vector & s, const SizeType & level) override
        // {
        //     if(level == 1)
        //     {
        //         _coarse_tr_subproblem->current_radius(get_delta(level-1));  
        //         _coarse_tr_subproblem->constrained_solve(H, g, s); 
        //     }
        //     else
        //     {
        //         _smoother_tr_subproblem->current_radius(get_delta(level-1));  
        //         _smoother_tr_subproblem->constrained_solve(H, g, s); 
        //     }

        //     return true; 
        // }





    };

}

#endif //UTOPIA_RMTR_INFTY_HPP

