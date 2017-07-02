/*
* @Author: alenakopanicakova
* @Date:   2017-04-19
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-02
*/

#ifndef UTOPIA_RMTR_HPP
#define UTOPIA_RMTR_HPP
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"

#include "utopia_TRSubproblem.hpp"
#include "utopia_Linear.hpp"
#include "utopia_Level.hpp"
#include "utopia_LS_Strategy.hpp"

#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_TR_base.hpp"

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
    class RMTR : public NonlinearMultiLevelBase<Matrix, Vector, FunctionType>,
                       public TrustRegionBase<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;
        typedef utopia::NonLinearSolver<Matrix, Vector>     Solver;
        typedef utopia::NonLinearSmoother<Matrix, Vector>   Smoother;
        typedef utopia::TRSubproblem<Matrix, Vector>        TRSubproblem; 
        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef utopia::Level<Matrix, Vector>               Level;

    
    public:

       /**
        * @brief      Multigrid class
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level. 
        */
        RMTR(    
                const std::shared_ptr<TRSubproblem> &tr_subproblem_coarse = std::shared_ptr<TRSubproblem>(),
                const std::shared_ptr<TRSubproblem> &tr_subproblem_smoother = std::shared_ptr<TRSubproblem>(),
                const Parameters params = Parameters()): 
                NonlinearMultiLevelBase<Matrix,Vector, FunctionType>(params), 
                _coarse_tr_subproblem(tr_subproblem_coarse), 
                _smoother_tr_subproblem(tr_subproblem_smoother) 
        {
            set_parameters(params); 
        }

        virtual ~RMTR(){} 
        

        void set_parameters(const Parameters params) override
        {
            NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::set_parameters(params);    
            _it_global                  = 0;          
            _delta_init                 = 1000000; 
            _parameters                 = params; 

            _max_coarse_it              = params.max_coarse_it();
            _max_smoothing_it           = params.max_smoothing_it();
            _eps_delta_termination      = params.eps_delta_termination();
            _delta_min                  = params.delta_min();
            _grad_smoothess_termination = params.grad_smoothess_termination();
            _eps_grad_termination       = params.eps_grad_termination();
            _hessian_update_delta       = params.hessian_update_delta();
            _hessian_update_eta         = params.hessian_update_eta();

        }

        using NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::solve; 

        virtual std::string name_id() override { return "RMTR";  }
        

        /**
         * @brief      The solve function for multigrid method. 
         *
         * @param[in]  rhs   The right hand side.
         * @param      x_0   The initial guess. 
         *
         */
        virtual bool solve(FunctionType &fine_fun, Vector & x_h, const Vector & rhs) override
        {
            
            Vector F_h  = local_zeros(local_size(x_h)); 

            bool converged = false; 
            SizeType l = this->num_levels(); 
            Scalar r_norm, r0_norm, rel_norm;
            
            std::cout<< this->name_id() <<"     Number of levels: "<< l << "  \n"; 

            Matrix hessian; 
            fine_fun.hessian(x_h, hessian); 


            fine_fun.gradient(x_h, F_h); 
            r0_norm = norm2(F_h); 

            this->make_iterate_feasible(fine_fun, x_h); 
            Scalar energy; 
            fine_fun.value(x_h, energy); 

            _it_global = 1; 


            //-------------- INITIALIZATIONS ---------------
            init_deltas(); 
            init_delta_gradients(); 
            init_x_initials(); 

            if(CONSISTENCY_LEVEL == SECOND_ORDER || CONSISTENCY_LEVEL == GALERKIN)
                init_delta_hessians(); 
            //----------------------------------------------


            while(!converged)
            {            
                if(this->cycle_type() =="multiplicative")
                    this->multiplicative_cycle(fine_fun, x_h, rhs, l); 
                else
                    std::cout<<"ERROR::UTOPIA_MG<< unknown MG type... \n"; 


                #ifdef CHECK_NUM_PRECISION_mode
                    if(has_nan_or_inf(x_h) == 1)
                    {
                        x_h = local_zeros(local_size(x_h));
                        return true; 
                    }
                #endif    

                fine_fun.gradient(x_h, F_h); 
                fine_fun.value(x_h, energy); 
                
                r_norm = norm2(F_h);
                rel_norm = r_norm/r0_norm; 

                ColorModifier red(FG_LIGHT_MAGENTA);
                ColorModifier def(FG_DEFAULT);
                std::cout << red; 
                std::cout<<"outer loop, it: "<< _it_global <<  "     ||g||:    "<< r_norm << "  energy: "<< energy << "   \n"; 
                std::cout << def; 

                // check convergence and print interation info
                converged = NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::check_convergence(_it_global, r_norm, rel_norm, 1); 
                _it_global++; 
            
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


        virtual bool multiplicative_cycle(FunctionType &fine_fun, Vector & u_l, const Vector &/*f*/, const SizeType & level) override
        {
            Vector g_fine, g_coarse, g_diff, g_restricted, u_2l, s_coarse, s_fine; 
            Matrix H_fine, H_restricted, H_coarse, H_diff; 

            Scalar ared=0.0, coarse_reduction=0.0, rho=0.0; 
            Scalar E_old, E_new; 
            bool converged = false; 


//--------------------------------------------------------------------------------------------------------------------------------------------
            // PRE-SMOOTHING 
            local_tr_solve(fine_fun, u_l, level); 
//--------------------------------------------------------------------------------------------------------------------------------------------

            fine_fun.gradient(u_l, g_fine);   

            if(level == this->num_levels())
            {
                converged =  this->criticality_measure_termination(norm2(g_fine)); 
                if(converged==true)
                    return true; 
            }

            // TODO:: check this out 
            // r_h = g_fine; // - f; 

            transfers(level-2).restrict(g_fine, g_restricted);
            transfers(level-2).project_down(u_l, u_2l); 

            this->make_iterate_feasible(levels(level-2), u_2l); 

            // here, grad should be zero by default on places where are BC conditions             
            levels(level-2).gradient(u_2l, g_coarse); 


            if(CONSISTENCY_LEVEL != GALERKIN)
                this->zero_boundary_correction(levels(level-2), g_restricted); 

            g_diff = g_restricted - g_coarse;  // tau correction 


            if(CONSISTENCY_LEVEL == SECOND_ORDER || CONSISTENCY_LEVEL == GALERKIN)
            {
                fine_fun.hessian(u_l, H_fine);   
                transfers(level-2).restrict(H_fine, H_restricted);
                
                if(CONSISTENCY_LEVEL == SECOND_ORDER)
                    this->zero_boundary_correction_mat(levels(level-2), H_restricted); 

                levels(level-2).hessian(u_2l, H_coarse); 
                H_diff = H_restricted - H_coarse; 
            }


            //----------------------------------------------------------------------------
            //                   initializing levels 
            //----------------------------------------------------------------------------
            set_delta(level-2, get_delta(level-1)); 
            set_delta_gradient(level-2, g_diff); 

            if(CONSISTENCY_LEVEL == SECOND_ORDER || CONSISTENCY_LEVEL == GALERKIN)
                set_delta_hessian(level-2, H_diff); 

            set_x_initial(level-2, u_2l); 
            set_delta_zero(level-2, get_delta(level-1)); 
        

            //----------------------------------------------------------------------------
            //               recursion  / Taylor correction
            //----------------------------------------------------------------------------

            // if grad is not smooth enoguh, we proceed to Taylor iterations, no recursion anymore
            if(level == 2 || grad_smoothess_termination(g_restricted, g_coarse))
            {
                SizeType l_new = level - 1; 
                coarse_reduction = local_tr_solve(levels(level-2), u_2l, l_new); 
            }
            else
            {
                // recursive call into RMTR 
                for(SizeType k = 0; k < this->mg_type(); k++)
                {   
                    SizeType l_new = level - 1; 
                                                        // check g_diff here !!!!!!!!!!!!!!!!!!!!!!!
                    multiplicative_cycle(levels(level-2), u_2l, g_diff, l_new); 
                }
            }

            //----------------------------------------------------------------------------
            //                       building trial point 
            //----------------------------------------------------------------------------

            s_coarse = u_2l - get_x_initial(level - 2);
            transfers(level-2).interpolate(s_coarse, s_fine);
            this->zero_boundary_correction(fine_fun, s_fine); 

            Vector u_t = u_l + s_fine; 

            Vector s_global = 0 * u_l; 
            if(level < this->num_levels())
                s_global = u_l - get_x_initial(level - 1);                  


            E_old = get_multilevel_energy(fine_fun,  u_l,  get_delta_gradient(level-1),  get_delta_hessian(level-1), s_global, level); 
            E_new = get_multilevel_energy(fine_fun,  u_t,  get_delta_gradient(level-1),  get_delta_hessian(level-1), s_global, level); 
            
            //----------------------------------------------------------------------------
            //                        trial point acceptance  
            //----------------------------------------------------------------------------

            ared = E_old - E_new; 
            rho = ared / coarse_reduction; 
            if(coarse_reduction<=0)
                rho = 0; 

            if(rho > this->rho_tol())
                u_l = u_t; 
            else
                std::cout<<"RMTR:: not taking trial point... \n"; 

            //----------------------------------------------------------------------------
            //                                  trust region update 
            //----------------------------------------------------------------------------
            delta_update(rho, level, s_global, converged); 
            if(converged==true && level == this->num_levels()){
                std::cout<<"YEs, second one ....... \n"; 
                return true; 
            }


            // just to see what is being printed 
            this->init_solver("RMTR_coarse_corr_stat", {" it. ", "   E_old     ", "   E_new", "ared   ",  "  coarse_level_reduction  ", "  rho  ", "  delta "}); 
            PrintInfo::print_iter_status(_it_global, {E_old, E_new, ared, coarse_reduction, rho, get_delta(level-1) }); 


    //--------------------------------------------------------------------------------------------------------------------------------------------
            // POST-SMOOTHING 
            local_tr_solve(fine_fun, u_l, level); 
    //--------------------------------------------------------------------------------------------------------------------------------------------


            return true; 
        }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        virtual Scalar local_tr_solve(FunctionType &fun, Vector & x, const SizeType & level)
        {   
            Vector g_diff, s_global, g; 
            Matrix H_diff, H; 

            if(level < this->num_levels())
            {
                g_diff = this->get_delta_gradient(level-1); 

                if(CONSISTENCY_LEVEL == SECOND_ORDER || CONSISTENCY_LEVEL == GALERKIN)
                    H_diff = this->get_delta_hessian(level-1); 
            }


            SizeType it_success = 0, it = 0; 
            Scalar ared = 0. , pred = 0., rho = 0., energy_old=9e9, energy_new=9e9, g_norm=1.0, reduction = 0.0; 
            bool make_grad_updates = true, make_hess_updates = true, converged = false; 


            Vector s = local_zeros(local_size(x)); 
            s_global = s; 


            get_multilevel_gradient(fun, x, g, g_diff, H_diff, s_global, level); 
            energy_old = get_multilevel_energy(fun,  x,  g_diff,  H_diff, s_global, level); 
            g_norm = norm2(g); 
            

            print_level_info(level); 
            PrintInfo::print_iter_status(0, {g_norm, energy_old, ared, pred, rho, get_delta(level-1) }); 

            it++; 

            while(!converged)
            {
                
                get_multilevel_hessian(fun, x, H, H_diff, level); 
                energy_old = get_multilevel_energy(fun,  x,  g_diff,  H_diff, s_global, level); 

            //----------------------------------------------------------------------------
            //     solving constrained system to get correction and  building trial point 
            //----------------------------------------------------------------------------
                // correction needs to get prepared 
                s = 0 * x;
                this->solve_qp_subproblem(H, g, s, level); 

                TrustRegionBase<Matrix, Vector>::get_pred(g, H, s, pred); 
                Vector tp = x + s;  
                
                if(level < this->num_levels())
                    s_global = tp - get_x_initial(level - 1);                  

                energy_new = get_multilevel_energy(fun,  tp,  g_diff,  H_diff, s_global, level); 
                ared = energy_old - energy_new; 

                // choice for the moment  - check for nans !!! 
                rho = ared/pred; 
                rho = (rho != rho) ? 0.0 : rho; 

            //----------------------------------------------------------------------------
            //     acceptance of trial point 
            //----------------------------------------------------------------------------
              
                // good reduction, accept trial point 
                if (rho >= this->rho_tol())
                {
                    x = tp; 
                    reduction += ared; 
                    it_success++; 
                    make_grad_updates =  true; 
                    // energy_old = energy_new; 
                }
                else
                {
                    // since point was not taken 
                    s_global -= s; 
                    make_grad_updates =  false; 
                }

            //----------------------------------------------------------------------------
            //     trust region update 
            //----------------------------------------------------------------------------
           
                this->delta_update(rho, level, s_global, converged); 

                if(make_grad_updates)
                {
                    Vector g_old = g; 
                    get_multilevel_gradient(fun, x, g, g_diff, H_diff, s_global, level); 
                    g_norm = norm2(g); 
                    make_hess_updates =  update_hessian(g, g_old, s, H, rho, g_norm); 
                }

                converged  = (converged  == true) ? true : check_local_convergence(it_success,  g_norm, level, get_delta(level-1)); 
                PrintInfo::print_iter_status(it, {g_norm, energy_new, ared, pred, rho, get_delta(level-1)}); 
                it++; 

            }

            ColorModifier color_def(FG_DEFAULT);
            std::cout<< color_def; 

            return reduction; 
        }





    protected:
        virtual void delta_update(const Scalar & rho, const SizeType & level, const Vector & s_global, bool & converged)
        {
            Scalar intermediate_delta; 

            if(rho < this->eta1())
            {   
                 intermediate_delta = this->gamma1() * get_delta(level-1); 
            }
            else if (rho > this->eta2() )
            {
                 intermediate_delta = this->gamma2() * get_delta(level-1); 
            }
            else
            {
                intermediate_delta = get_delta(level-1); 
            }      


            // on the finest level we work just with one radius 
            if(level==this->num_levels())
            {
                set_delta(level-1, intermediate_delta); 
            }
            else
            {
                Scalar corr_norm = level_dependent_norm(s_global, level); 
                converged = delta_termination(corr_norm, level); 
                if(converged)
                {
                    std::cout<<"termination  due to small radius for given level ... \n"; 
                    return; 
                }

                corr_norm = get_delta_zero(level-1) - corr_norm; 
                corr_norm = std::min(intermediate_delta, corr_norm); 

                if(corr_norm <= 0)
                    corr_norm = 0; 

                set_delta(level-1, std::min(intermediate_delta, corr_norm)); 
            }
        }



        virtual bool check_local_convergence(const SizeType & it_success, const Scalar & g_norm, const SizeType & level, const Scalar & delta)
        {   
            // coarse one 
            if(level == 1 && (it_success >= _max_coarse_it))
                return true; 
            // every other level 
            else if (level > 1 && it_success >= _max_smoothing_it)
                return true; 

            if(delta < _delta_min)
                return true; 

            return criticality_measure_termination(g_norm); 
        }

        /**
         * @brief      THis check guarantee that iterates at a lower level remain in the TR radius defined at the calling level
         *
         * @param[in]  corr_norm  The norm of sum of all corrections on given level 
         * @param[in]  level      The level
         *
         */
        virtual bool delta_termination(const Scalar & corr_norm, const SizeType & level)
        {   
            return (corr_norm > (1 - _eps_delta_termination) * get_delta(level-1)) ? true : false; 
        }



        // needs to be overriden in case of \infty norm 
        virtual bool criticality_measure_termination(const Scalar & g_norm)
        {
            // this should be fancier based on Graffon paper, but it is quite boring to work with so many mesh informations
            // Scalar _eps_grad_termination = std::min(0.001, eps_grad_termination[level]/ psi); 
            
            return (g_norm < _eps_grad_termination) ? true : false;    
        }



        virtual bool grad_smoothess_termination(const Vector & g_restricted, const Vector & g_coarse)
        {
            Scalar Rg_norm = norm2(g_restricted); 
            Scalar g_norm = norm2(g_coarse);
            // return (Rg_norm >= _grad_smoothess_termination g_norm) ? true : false;   
            
            if(Rg_norm >= _grad_smoothess_termination * g_norm)
            {
                std::cout<<"grad is not smooth enough............... => NO RECURSION ANYMORE ....  \n"; 
                return true; 
            }
            else
                return false; 

        }



        virtual bool update_hessian(const Vector & g_new, const Vector & g_old, const Vector & s, const Matrix & H, const Scalar & rho, const Scalar & g_norm)
        {
            // iteration is not sucessful enough
            if(rho > 0 && rho < _hessian_update_eta)
                return true; 

            Vector help = g_new - g_old - H * s; 

            // Hessian approx is relativelly poor
            return (norm2(help) > _hessian_update_delta * g_norm) ? true : false; 
        }




        virtual bool set_delta(const SizeType & level, const Scalar & radius)
        {
            _deltas[level] = radius; 
            return true; 
        }


        virtual Scalar get_delta(const SizeType & level) const 
        {
            return _deltas[level]; 
        }



        virtual bool set_delta_zero(const SizeType & level, const Scalar & radius)
        {
            _deltas_zero[level] = radius; 
            return true; 
        }


        virtual Scalar get_delta_zero(const SizeType & level) const 
        {
            return _deltas_zero[level]; 
        }



        // organized from coarsest
        // delta[0] =  coarsest level
        virtual bool init_deltas()
        {
            for(Scalar i = 0; i < this->num_levels(); i ++)
                _deltas.push_back(_delta_init); 

            for(Scalar i = 0; i < this->num_levels()-1; i ++)
                _deltas_zero.push_back(0); 

            return true; 
        }


        // organized from coarsest
        // _delta_hessians[0] =  coarsest level
        // NOTE: we do not have any for the finest level, since function on the finest level is taken from problem definition by itself
        bool init_delta_hessians()
        {
            _delta_hessians.resize(this->num_levels()-1); 
            return true; 
        }


        virtual bool set_delta_hessian(const SizeType & level, const Matrix & H_diff)
        {
            _delta_hessians[level] = H_diff; 
            return true; 
        }


        virtual Matrix & get_delta_hessian(const SizeType & level) 
        {
            return _delta_hessians[level]; 
        }


        // organized from coarsest
        // _x_initials[0] =  coarsest level
        // NOTE: we do not have any for the finest level, since function on the finest level is taken from problem definition by itself
        virtual bool init_x_initials()
        {
            _x_initials.resize(this->num_levels()-1); 
            return true; 
        }


        virtual bool set_x_initial(const SizeType & level, const Vector & x)
        {
            _x_initials[level] = x; 
            return true; 
        }


        virtual Vector & get_x_initial(const SizeType & level) 
        {
            return _x_initials[level]; 
        }


        // organized from coarsest
        // _delta_gradients[0] =  coarsest level
        // NOTE: we do not have any for the finest level, since function on the finest level is taken from problem definition by itself
        virtual bool init_delta_gradients()
        {
            _delta_gradients.resize(this->num_levels()-1); 
            return true; 
        }


        virtual bool set_delta_gradient(const SizeType & level, const Vector & g_diff)
        {
            _delta_gradients[level] = g_diff; 
            return true; 
        }


        virtual Vector & get_delta_gradient(const SizeType & level) 
        {
            return _delta_gradients[level]; 
        }



        // in order to be able to use full cycle 
        // CHECK IF RHS does not need to be set-up
        virtual bool coarse_solve(FunctionType &fun, Vector &x, const Vector & /*rhs*/) override
        {
            local_tr_solve(fun, x, 0); 
            return true; 
        }



        virtual void print_level_info(const SizeType & level)
        {
            ColorModifier color_out(FG_LIGHT_YELLOW);
            if(this->verbose())
            {
                if(level == 1)
                {
                    std::cout << color_out; 
                    this->init_solver("COARSE SOLVE", {" it. ", "|| g_norm ||", "   E + <g_diff, s>", "ared   ",  "  pred  ", "  rho  ", "  delta "}); 
                }
                else
                {
                    color_out.set_color_code(FG_LIGHT_GREEN); 
                    std::cout << color_out; 
                    this->init_solver("SMOOTHER", {" it. ", "|| g_norm ||", "   E + <g_diff, s>", "ared   ",  "  pred  ", "  rho  ", "  delta "}); 
                }

                std::cout<<"LEVEL: "<< level << "   \n"; 
            }
        }



        virtual bool solve_qp_subproblem(const Matrix & H, const Vector & g, Vector & s, const SizeType & level)
        {
            if(level == 1)
            {
                _coarse_tr_subproblem->current_radius(get_delta(level-1));  
                _coarse_tr_subproblem->tr_constrained_solve(H, g, s); 
            }
            else
            {
                _smoother_tr_subproblem->current_radius(get_delta(level-1));  
                _smoother_tr_subproblem->tr_constrained_solve(H, g, s); 
            }

            return true; 
        }



        virtual bool get_multilevel_hessian(const FunctionType & fun, const Vector & x,  Matrix & H, const Matrix & H_diff, const SizeType & level)
        {
            if(level < this->num_levels())
                return MultilevelHessianEval<Matrix, Vector, FunctionType, CONSISTENCY_LEVEL>::compute_hessian(fun, x, H, H_diff);
            else
                return fun.hessian(x, H); 
        }


        virtual bool get_multilevel_gradient(const FunctionType & fun, const Vector & x,  Vector & g, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global, const SizeType & level)
        {
            if(level < this->num_levels())
                return MultilevelGradientEval<Matrix, Vector, FunctionType, CONSISTENCY_LEVEL>::compute_gradient(fun, x, g, g_diff, H_diff, s_global);
            else
                 return fun.gradient(x, g); 
        }


        virtual Scalar get_multilevel_energy(const FunctionType & fun, const Vector & x, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global, const SizeType & level)
        {
            if(level < this->num_levels())
                return MultilevelEnergyEval<Matrix, Vector, FunctionType, CONSISTENCY_LEVEL>::compute_energy(fun, x, g_diff, H_diff, s_global); 
            else
            {
                Scalar energy; 
                fun.value(x, energy); 
                return energy; 
            }
        }



private: 
        virtual  Scalar level_dependent_norm(const Vector & u, const SizeType & current_l)
        {
            if(current_l == this->num_levels())
                return 0.0; 
            else
            {
                Vector s = u; // carries over prolongated correction
         
                for(SizeType i = current_l; i < this->num_levels(); i++)
                    transfers(i-1).interpolate(s, s); 
                return norm2(s); 
            }    
        }



    protected:   
        SizeType                            _it_global; 
        std::vector<Scalar>                 _deltas;  
        std::vector<Scalar>                 _deltas_zero;        

        


        std::shared_ptr<TRSubproblem>        _coarse_tr_subproblem; 
        std::shared_ptr<TRSubproblem>        _smoother_tr_subproblem; 


        std::vector<Vector>           _delta_gradients; 
        std::vector<Matrix>           _delta_hessians; 
        std::vector<Vector>           _x_initials; 


        // ----------------------- PARAMETERS ----------------------
        Parameters                          _parameters; 
        Scalar                              _delta_init; 
        SizeType                            _max_coarse_it; 
        SizeType                            _max_smoothing_it; 

        Scalar                         _eps_delta_termination; 
        Scalar                         _delta_min; 

        Scalar                         _grad_smoothess_termination; 
        Scalar                         _eps_grad_termination; 

        Scalar                         _hessian_update_delta; 
        Scalar                         _hessian_update_eta; 


    };

}

#endif //UTOPIA_RMTR_HPP

