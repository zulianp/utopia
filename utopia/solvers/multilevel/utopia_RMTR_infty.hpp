/*
* @Author: alenakopanicakova
* @Date:   2017-04-19
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-06-10
*/

#ifndef UTOPIA_RMTR_INFTY_HPP
#define UTOPIA_RMTR_INFTY_HPP
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"

#include "utopia_TRSubproblem.hpp"
#include "utopia_Linear.hpp"
#include "utopia_Level.hpp"


#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"


namespace utopia 
{
    /**
     * @brief      The class for Nonlinear Multigrid solver. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector, class FunctionType>
    class RMTR_infty : public NonlinearMultiLevelBase<Matrix, Vector, FunctionType>,
                       public TrustRegionBase<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::NonLinearSolver<Matrix, Vector>     Solver;
        typedef utopia::NonLinearSmoother<Matrix, Vector>   Smoother;
        typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblem; 
        typedef utopia::Transfer<Matrix, Vector>   Transfer;


        typedef utopia::Level<Matrix, Vector>               Level;

    

        enum Coherence_level {  FIRST_ORDER  = 1, 
                                SECOND_ORDER = 2, 
                                GALERKIN     = 0};


    public:

       /**
        * @brief      Multigrid class
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level. 
        */
        RMTR_infty(    
                const std::shared_ptr<TRSubproblem> &tr_subproblem_coarse = std::shared_ptr<TRSubproblem>(),
                const std::shared_ptr<TRSubproblem> &tr_subproblem_smoother = std::shared_ptr<TRSubproblem>(),
                const Parameters params = Parameters()): 
                NonlinearMultiLevelBase<Matrix,Vector, FunctionType>(params), 
                _coarse_tr_subproblem(tr_subproblem_coarse), 
                _smoother_tr_subproblem(tr_subproblem_smoother)
        {
            set_parameters(params); 
        }

        virtual ~RMTR_infty(){} 
        

        void set_parameters(const Parameters params)  // override
        {
            NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::set_parameters(params);             
            _delta_init     = 1000000; 
            _parameters     = params; 
            _coherence      = FIRST_ORDER, 
            _max_coarse_it  = 10;  
            _max_fine_it    = 2;
        }

        virtual bool solve(FunctionType & fine_fun, Vector &x_h)
        {
            Vector rhs = local_zeros(local_size(x_h)); 
            return solve_rhs(fine_fun,  x_h, rhs); 
        }

        /**
         * @brief      The solve function for multigrid method. 
         *
         * @param[in]  rhs   The right hand side.
         * @param      x_0   The initial guess. 
         *
         */
        virtual bool solve_rhs(FunctionType &fine_fun, Vector & x_h, const Vector & rhs) 
        {
            
            Vector F_h  = local_zeros(local_size(x_h)); 

            bool converged = false; 
            SizeType it = 0, l = this->num_levels(); 
            Scalar r_norm, r0_norm, rel_norm;
            std::cout<<"RMTR_infty: number of levels: "<< l << "  \n"; 

            Matrix hessian; 
            fine_fun.hessian(x_h, hessian); 


            fine_fun.gradient(x_h, F_h); 
            r0_norm = norm2(F_h); 

            this->make_iterate_feasible(fine_fun, x_h); 
            Scalar energy; 
            fine_fun.value(x_h, energy); 

            it++; 


            //-------------- INITIALIZATIONS ---------------
            // init deltas 
            init_deltas(); 
            init_delta_gradients(); 
            init_x_initials(); 

            if(_coherence == SECOND_ORDER || _coherence == GALERKIN)
                init_delta_hessians(); 
            //----------------------------------------------


            while(!converged)
            {            
                if(this->cycle_type() =="multiplicative")
                    multiplicative_cycle(fine_fun, x_h, rhs, l, it); 
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
                std::cout<<"outer loop, it: "<< it <<  "     ||g||:    "<< r_norm << "  energy: "<< energy << "   \n"; 
                std::cout << def; 

                // check convergence and print interation info
                converged = NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::check_convergence(it, r_norm, rel_norm, 1); 
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


        bool multiplicative_cycle(FunctionType &fine_fun, Vector & u_l, const Vector &f, const SizeType & level, const SizeType & it_global)
        {
            Vector g_fine, g_coarse, g_diff, r_h,  g_restricted, u_2l, s_coarse, s_fine, u_init; 
            Matrix H_fine, H_restricted, H_coarse, H_diff; 

            // tr radius on this level 
            Scalar delta = get_delta(level-2); 


            Scalar ared=0.0, coarse_reduction=0.0, rho=0.0; 

//--------------------------------------------------------------------------------------------------------------------------------------------
            // PRE-SMOOTHING 
            local_tr_solve(fine_fun, u_l,  s_coarse, coarse_reduction, level); 
//--------------------------------------------------------------------------------------------------------------------------------------------

            fine_fun.gradient(u_l, g_fine);   

            r_h = g_fine - f; 

            transfers(level-2).restrict(r_h, g_restricted);
            transfers(level-2).project_down(u_l, u_2l); 

            this->make_iterate_feasible(levels(0), u_2l); 

            // here, grad should be zero by default on places where are BC conditions             
            levels(level-2).gradient(u_2l, g_coarse); 

            u_init = u_2l; 

            if(_coherence != GALERKIN)
                this->zero_boundary_correction(levels(0), g_restricted); 

            g_diff = g_restricted - g_coarse;  // tau correction 


            if(_coherence == SECOND_ORDER || _coherence == GALERKIN)
            {
                fine_fun.hessian(u_l, H_fine);   
                transfers(level-2).restrict(H_fine, H_restricted);
                
                if(_coherence == SECOND_ORDER)
                    this->zero_boundary_correction_mat(levels(0), H_restricted); 

                levels(level-2).hessian(u_2l, H_coarse); 
                H_diff = H_restricted - H_coarse; 
            }


            if(level == 2)
            {
                SizeType l_new = level - 1; 
                set_delta(0, get_delta(l_new)); 

                set_delta_gradient(l_new - 1, g_diff); 

                if(_coherence == SECOND_ORDER || _coherence == GALERKIN)
                    set_delta_hessian(l_new - 1, H_diff); 

                set_x_initial(l_new - 1, u_2l); 
                local_tr_solve(levels(0), u_2l, s_coarse, coarse_reduction, l_new); 
            }
            else
            {
                // // recursive call into FAS - needs to be checked 
                // for(SizeType k = 0; k < this->mg_type(); k++)
                // {   
                //     SizeType l_new = l - 1; 
                //     multiplicative_cycle(levels(l-2), u_2l, L_2l, l_new); 
                // }
            }


            transfers(level-2).interpolate(s_coarse, s_fine);
            this->zero_boundary_correction(fine_fun, s_fine); 

            Scalar E_old, E_new; 
            fine_fun.value(u_l, E_old); 

            Vector u_t = u_l; 

            u_t += s_fine; 
            fine_fun.value(u_t, E_new); 

            ared = E_old - E_new; 
            rho = ared / coarse_reduction; 

           std::cout<<"E_old:  "<< E_old << "  E_new:  "<< E_new << "  \n"; 

            // TODO:: deltas, line-search, 

            if(coarse_reduction<=0)
                rho = 0; 

            if(rho > this->rho_tol())
                u_l = u_t; 
            else{
                std::cout<<"RMTR:: not taking trial point... \n"; 
            }


            //----------------------------------------------------------------------------
            //     trust region update 
            //----------------------------------------------------------------------------
           
           // TODO:: fix this ... for more levels 
            Vector s_global = 0 * u_l; 
            Scalar delta0 = 0.0; 
            delta_update(rho, level, delta0, s_global); 



            // just to see what is being printed 
            this->init_solver("RMTR_coarse_corr_stat", {" it. ", "   E_old     ", "   E_old", "ared   ",  "  coarse_reduction  ", "  rho  ", "  delta "}); 
            PrintInfo::print_iter_status(it_global, {E_old, E_new, ared, coarse_reduction, rho, get_delta(level-1) }); 




    //--------------------------------------------------------------------------------------------------------------------------------------------
            // POST-SMOOTHING 
            local_tr_solve(fine_fun, u_l, s_coarse, coarse_reduction, level); 
    //--------------------------------------------------------------------------------------------------------------------------------------------




            return true; 
        }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        bool local_tr_solve(FunctionType &fun, Vector & x,  Vector & s_global, Scalar & reduction, const SizeType & level)
        {   
            Vector g_diff; 
            Matrix H_diff; 

            if(level < this->num_levels())
            {
                g_diff = get_delta_gradient(level-1); 
                
                if(_coherence == SECOND_ORDER || _coherence == GALERKIN)
                {
                     H_diff = get_delta_hessian(level-1); 
                }
            }


            if(_coherence == FIRST_ORDER)
                std::cout<<"-------- Yes, first order....... \n"; 
            else if(_coherence == SECOND_ORDER)
                std::cout<<"-------- Yes, second order....... \n"; 

            std::cout<<"level: "<< level << "   \n"; 

            SizeType it_success = 0, it = 0; 
            bool converged = false; 
            bool make_updates = true; 


            Scalar energy_old, energy_new, energy_init, g_norm;
            Scalar ared = 0 , pred = 0, rho = 0; 

            ColorModifier color_out(FG_LIGHT_YELLOW);
            ColorModifier color_def(FG_DEFAULT);

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

            Scalar delta0 = get_delta(level-1); 

            std::cout<<"coarse delta0: "<< delta0 << "  \n"; 


            Vector s = 0 * x; 
            s_global = s; 
            reduction = 0.0; 


            Vector g; 
            Matrix H; 

            
            // --------------------------------------- computation of grad -------------------------------
            if(_coherence != GALERKIN)
                fun.gradient(x, g);

            if(level < this->num_levels())
            {
                if(_coherence == FIRST_ORDER)
                    g += g_diff; 
                else if(_coherence == SECOND_ORDER)
                    g += g_diff + H_diff * s_global; 
                else if(_coherence == GALERKIN)
                    g = g_diff + H_diff * s_global; 
            }
            // --------------------------------------------------------------------------------------------


            g_norm = norm2(g); 

            // --------------------------------------- computation of energy -------------------------------
            if(_coherence != GALERKIN)
                fun.value(x, energy_init); 

            if(level < this->num_levels())
            {
                if(_coherence == FIRST_ORDER)
                    energy_init += dot(g_diff, s_global); 
                else if(_coherence == SECOND_ORDER)
                    energy_init += dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                else if(_coherence == GALERKIN)
                    energy_init = dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
            }
            // --------------------------------------------------------------------------------------------
            


            PrintInfo::print_iter_status(0, {g_norm, energy_init, ared, pred, rho, get_delta(level-1) }); 

            it++; 

            while(!converged)
            {
                
                // --------------------------------------- computation of hessian -------------------------------
                if(make_updates)
                {
                    if(_coherence != GALERKIN)
                        fun.hessian(x, H); 

                    if(level < this->num_levels())
                    {
                        if(_coherence == SECOND_ORDER)
                            H = H + H_diff; 
                        else if(_coherence == GALERKIN)
                            H = H_diff; 
                    }
                }
                // --------------------------------------------------------------------------------------------

                // --------------------------------------- computation of energy -------------------------------
                if(make_updates)
                {
                    if(_coherence != GALERKIN)
                        fun.value(x, energy_old); 

                    if(level < this->num_levels())
                    {
                        if(_coherence == FIRST_ORDER)
                            energy_old += dot(g_diff, s_global); 
                        else if(_coherence == SECOND_ORDER)
                            energy_old += dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                        else if(_coherence == GALERKIN)
                            energy_old = dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                    }
                }
                // --------------------------------------------------------------------------------------------


            //----------------------------------------------------------------------------
            //     solving constrained system to get correction
            //----------------------------------------------------------------------------
                // this needs to get prepared 
                s = 0 * x;

                if(level == 1)
                {
                    _coarse_tr_subproblem->current_radius(get_delta(level-1));  
                    _coarse_tr_subproblem->constrained_solve(H, g, s); 
                }
                else
                {
                    _smoother_tr_subproblem->current_radius(get_delta(level-1));  
                    _smoother_tr_subproblem->constrained_solve(H, g, s); 
                }

                TrustRegionBase<Matrix, Vector>::get_pred(g, H, s, pred); 

            //----------------------------------------------------------------------------
            //     building trial point 
            //----------------------------------------------------------------------------
                
                Vector tp = x + s;  
                // s_global += s;  
                
                if(level < this->num_levels())
                {
                    s_global = tp - get_x_initial(level - 1);                  
                }

            // --------------------------------------- computation of energy -------------------------------
                if(_coherence != GALERKIN)
                    fun.value(tp, energy_new); 

                if(level < this->num_levels())
                {
                    if(_coherence == FIRST_ORDER)
                        energy_new += dot(g_diff, s_global); 
                    else if(_coherence == SECOND_ORDER)
                        energy_new += dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                    else if(_coherence == GALERKIN)
                        energy_new = dot(g_diff, s_global) + 0.5 * dot(s_global, H_diff * s_global); 
                }
            // --------------------------------------------------------------------------------------------


                ared = energy_old - energy_new; 

                // choice for the moment 
                rho = ared/pred; 

                //----------------------------------------------------------------------------
                //     acceptance of trial point 
                //----------------------------------------------------------------------------
                  
                  // good reduction, accept trial point 
                  if (rho >= this->rho_tol())
                  {
                    x = tp; 
                    reduction += ared; 
                    it_success++; 
                    make_updates =  true; 
                  }
                  else
                  {
                    // since point was not taken 
                    s_global -= s; 
                    make_updates =  false; 
                  }

                //----------------------------------------------------------------------------
                //     trust region update 
                //----------------------------------------------------------------------------
               
                delta_update(rho, level, delta0, s_global); 


                // --------------------------------------- computation of grad -------------------------------
                if(make_updates)
                {
                    if(_coherence != GALERKIN)
                        fun.gradient(x, g);

                    if(level < this->num_levels())
                    {
     
                        if(_coherence == FIRST_ORDER)
                            g += g_diff; 
                        else if(_coherence == SECOND_ORDER)
                            g += g_diff + H_diff * s_global; 
                        else if(_coherence == GALERKIN)
                            g = g_diff + H_diff * s_global; 
                    }

                    g_norm = norm2(g); 
                }
                // --------------------------------------------------------------------------------------------

                converged  = check_convergence(it_success,  g_norm, level, get_delta(level-1)); 
                PrintInfo::print_iter_status(it, {g_norm, energy_new, ared, pred, rho, get_delta(level-1)}); 
                it++; 

            }

            std::cout<< color_def; 

            return true; 
        }





    public:
        Scalar level_dependent_norm(const Vector & u, const SizeType & current_l)
        {
            if(current_l == this->num_levels())
            {
                return 0.0; 
            }
            else
            {
                Vector s = u; // carries over prolongated correction
         
                for(SizeType i = current_l; i < this->num_levels(); i++)
                {   
                    transfers(i-1).interpolate(s, s); 
                }
                return norm2(s); 
            }    
        }


        virtual void delta_update(const Scalar & rho, const SizeType & level, const Scalar & delta0, const Vector & s_global)
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
                std::cout<<"bigger delta   \n"; 
                set_delta(level-1, intermediate_delta); 
            }
            else
            {
                std::cout<<"smaller delta ... \n"; 
                // compute radius 
                Scalar delta_help = delta0 - level_dependent_norm(s_global, level); 
                // delta = std::min(intermediate_delta, delta_help); 
                set_delta(level-1, std::min(intermediate_delta, delta_help)); 
            }
        }



        bool check_convergence(const SizeType & it_success, const Scalar & g_norm, const SizeType & level, const Scalar & delta)
        {   

            // coarse one 
            if(level == 1)
            {
                if(it_success >= _max_coarse_it)
                    return true; 
            }
            // every other level 
            else
            {
                if(it_success >= _max_fine_it)
                    return true; 
            }

            if(delta < 1e-10)
                return true; 


            if(g_norm <= 1e-8)
                return true; 
            else
                return false; 

        }

        // organized from coarsest
        // _delta_gradients[0] =  coarsest level
        // NOTE: we do not have any for the finest level, since function on the finest level is taken from problem definition by itself
        bool init_delta_gradients()
        {
            _delta_gradients.resize(this->num_levels()-1); 
            return true; 
        }


        bool set_delta_gradient(const SizeType & level, const Vector & g_diff)
        {
            _delta_gradients[level] = g_diff; 
            return true; 
        }


        Vector & get_delta_gradient(const SizeType & level) 
        {
            return _delta_gradients[level]; 
        }


        bool set_delta(const SizeType & level, const Scalar & radius)
        {
            _delta[level] = radius; 
            return true; 
        }


        Scalar get_delta(const SizeType & level) const 
        {
            return _delta[level]; 
        }


        // organized from coarsest
        // delta[0] =  coarsest level
        bool init_deltas()
        {
            for(Scalar i = 0; i < this->num_levels(); i ++)
                _delta.push_back(_delta_init); 

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


        bool set_delta_hessian(const SizeType & level, const Matrix & H_diff)
        {
            _delta_hessians[level] = H_diff; 
            return true; 
        }


        Matrix & get_delta_hessian(const SizeType & level) 
        {
            return _delta_hessians[level]; 
        }




        // organized from coarsest
        // _x_initials[0] =  coarsest level
        // NOTE: we do not have any for the finest level, since function on the finest level is taken from problem definition by itself
        bool init_x_initials()
        {
            _x_initials.resize(this->num_levels()-1); 
            return true; 
        }


        bool set_x_initial(const SizeType & level, const Vector & x)
        {
            _x_initials[level] = x; 
            return true; 
        }


        Vector & get_x_initial(const SizeType & level) 
        {
            return _x_initials[level]; 
        }





    protected:   
        Scalar                              _delta_init; 
        std::vector<Scalar>                 _delta;                  
        Coherence_level                     _coherence; 

    private:
        Parameters                          _parameters; 


        std::shared_ptr<TRSubproblem>        _coarse_tr_subproblem; 
        std::shared_ptr<TRSubproblem>        _smoother_tr_subproblem; 

        SizeType                            _max_coarse_it; 
        SizeType                            _max_fine_it; 


        std::vector<Vector>           _delta_gradients; 
        std::vector<Matrix>           _delta_hessians; 

        std::vector<Vector>           _x_initials; 


    };

}

#endif //UTOPIA_RMTR_INFTY_HPP

