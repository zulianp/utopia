/*
* @Author: alenakopanicakova
* @Date:   2017-04-19
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2018-02-08
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
#include "utopia_TRBase.hpp"

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
            _parameters                 = params; 

            _delta_init                 = params.delta0(); 
            _max_coarse_it              = params.max_coarse_it();
            _max_smoothing_it           = params.max_smoothing_it();
            _eps_delta_termination      = params.eps_delta_termination();
            _delta_min                  = params.delta_min();
            _grad_smoothess_termination = params.grad_smoothess_termination();
            _eps_grad_termination       = params.eps_grad_termination();
            _hessian_update_delta       = params.hessian_update_delta();
            _hessian_update_eta         = params.hessian_update_eta();

            _verbosity_level           = params.verbosity_level(); 

        }


        VerbosityLevel verbosity_level() const 
        {
            return _verbosity_level; 
        }

        void verbosity_level(const VerbosityLevel & level )
        {
            _verbosity_level = level; 
        }


        using NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::solve; 

        virtual std::string name_id() override { return "RMTR";  }
        

        void set_eps_grad_termination(const Scalar & eps_grad_termination)
        {
            _eps_grad_termination = eps_grad_termination; 
        }


        void set_max_coarse_it(const SizeType & max_coarse_it)
        {
            _max_coarse_it = max_coarse_it; 
        }


        void set_max_smoothing_it(const SizeType & max_smoothing_it)
        {
            _max_smoothing_it = max_smoothing_it; 
        }


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
            
            if(this->verbose())
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

            if(this->verbose())
            {
                ColorModifier red(FG_LIGHT_MAGENTA);
                ColorModifier def(FG_DEFAULT);
                std::cout << red; 

                if(verbosity_level() >= VERBOSITY_LEVEL_NORMAL)
                    this->init_solver("RMTR OUTER SOLVE", {" it. ", "|| g_norm ||", "   E "}); 

                PrintInfo::print_iter_status(_it_global, {r0_norm, energy}); 
                std::cout << def; 
            }


            while(!converged)
            {            
                if(this->cycle_type() == MULTIPLICATIVE_CYCLE)
                    this->multiplicative_cycle(fine_fun, x_h, rhs, l); 
                else{
                    std::cout<<"ERROR::UTOPIA_RMTR << unknown cycle type, solving in multiplicative manner ... \n"; 
                    this->multiplicative_cycle(fine_fun, x_h, rhs, l); 
                }


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

                if(this->verbose())
                {
                    ColorModifier red(FG_LIGHT_MAGENTA);
                    ColorModifier def(FG_DEFAULT);
                    std::cout << red; 

                    if(verbosity_level() > VERBOSITY_LEVEL_NORMAL)
                        this->init_solver("RMTR OUTER SOLVE", {" it. ", "|| g_norm ||", "   E "}); 

                    PrintInfo::print_iter_status(_it_global, {r_norm, energy}); 
                    std::cout << def; 
                }

                // check convergence and print interation info
                converged = NonlinearMultiLevelBase<Matrix, Vector, FunctionType>::check_convergence(_it_global, r_norm, rel_norm, 1); 
                converged = (converged==true || this->get_delta(l-1) < 1e-14) ? true : false; 

                _it_global++; 
            
            }

            // benchmarking
            print_statistics(); 


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


        /**
         * @brief      Multiplicative cycle 
         *
         * @param      fine_fun   The fine fun
         * @param      u_l        Current iterate
         * @param[in]  f          The right hand side
         * @param[in]  level      The level
         *
         */
        virtual bool multiplicative_cycle(FunctionType &fine_fun, Vector & u_l, const Vector &/*f*/, const SizeType & level) override
        {
            Vector g_fine, g_coarse, g_diff, g_restricted, u_2l, s_coarse, s_fine; 
            Matrix H_fine, H_restricted, H_coarse, H_diff; 

            Scalar ared=0.0, coarse_reduction=0.0, rho=0.0; 
            Scalar E_old, E_new; 
            bool converged = false; 


//--------------------------------------------------------------------------------------------------------------------------------------------
            // PRE-SMOOTHING 
            this->local_tr_solve(fine_fun, u_l, level); 
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
                this->zero_correction_related_to_equality_constrain(levels(level-2), g_restricted); 

            g_diff = g_restricted - g_coarse;  // tau correction 


            if(CONSISTENCY_LEVEL == SECOND_ORDER || CONSISTENCY_LEVEL == GALERKIN)
            {
                fine_fun.hessian(u_l, H_fine);   
                transfers(level-2).restrict(H_fine, H_restricted);
                
                if(CONSISTENCY_LEVEL == SECOND_ORDER)
                    this->zero_correction_related_to_equality_constrain_mat(levels(level-2), H_restricted); 

                levels(level-2).hessian(u_2l, H_coarse); 
                H_diff = H_restricted - H_coarse; 
            }


            //----------------------------------------------------------------------------
            //                   initializing levels 
            //----------------------------------------------------------------------------
            this->set_delta(level-2, get_delta(level-1)); 
            this->set_delta_gradient(level-2, g_diff); 

            if(CONSISTENCY_LEVEL == SECOND_ORDER || CONSISTENCY_LEVEL == GALERKIN)
                this->set_delta_hessian(level-2, H_diff); 

            this->set_x_initial(level-2, u_2l); 
            this->set_delta_zero(level-2, get_delta(level-1)); 
        

            //----------------------------------------------------------------------------
            //               recursion  / Taylor correction
            //----------------------------------------------------------------------------

            // TODO:: is this correct??? 
            // if grad is not smooth enoguh, we proceed to Taylor iterations, no recursion anymore
            // if(level == 2 || this->grad_smoothess_termination(g_restricted, g_fine))
            if(level == 2)
            {
                SizeType l_new = level - 1; 
                coarse_reduction = this->local_tr_solve(levels(level-2), u_2l, l_new); 
            }
            else
            {
                // recursive call into RMTR 
                for(SizeType k = 0; k < this->mg_type(); k++)
                {   
                    SizeType l_new = level - 1; 
                    
                    // coarse_reduction is missing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                    this->multiplicative_cycle(levels(level-2), u_2l, g_diff, l_new); 
                }
            }

            //----------------------------------------------------------------------------
            //                       building trial point 
            //----------------------------------------------------------------------------

            s_coarse = u_2l - this->get_x_initial(level - 2);
            transfers(level-2).interpolate(s_coarse, s_fine);
            this->zero_correction_related_to_equality_constrain(fine_fun, s_fine); 

            Vector u_t = u_l + s_fine; 

            Vector s_global = 0 * u_l; 
            if(level < this->num_levels())
                s_global = u_l - this->get_x_initial(level - 1);                  


            E_old = this->get_multilevel_energy(fine_fun,  u_l,  get_delta_gradient(level-1),  get_delta_hessian(level-1), s_global, level); 
            E_new = this->get_multilevel_energy(fine_fun,  u_t,  get_delta_gradient(level-1),  get_delta_hessian(level-1), s_global, level); 
            
            //----------------------------------------------------------------------------
            //                        trial point acceptance  
            //----------------------------------------------------------------------------

            ared = E_old - E_new; 
            rho = ared / coarse_reduction; 
            if(coarse_reduction<=0)
                rho = 0; 

            if(rho > this->rho_tol())
                u_l = u_t; 
            else{
                if(this->verbose() &&  verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE)
                    std::cout<<"RMTR:: not taking trial point... \n"; 
            }

            //----------------------------------------------------------------------------
            //                                  trust region update 
            //----------------------------------------------------------------------------
            this->delta_update(rho, level, s_global, converged); 
            if(converged==true && level == this->num_levels()){
                if(this->verbose() &&  verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE)
                    std::cout<<"YEs, second one ....... \n"; 
                return true; 
            }

            if(this->verbose() &&  verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE)
            {
                // just to see what is being printed 
                this->init_solver("RMTR_coarse_corr_stat", {" it. ", "   E_old     ", "   E_new", "ared   ",  "  coarse_level_reduction  ", "  rho  ", "  delta "}); 
                PrintInfo::print_iter_status(_it_global, {E_old, E_new, ared, coarse_reduction, rho, get_delta(level-1) }); 
            }


    //--------------------------------------------------------------------------------------------------------------------------------------------
            // POST-SMOOTHING 
            this->local_tr_solve(fine_fun, u_l, level); 
    //--------------------------------------------------------------------------------------------------------------------------------------------


            return true; 
        }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        /**
         * @brief      Local TR solve 
         *
         * @param      fun    Evaluation fnction
         * @param      x      Current iterate 
         * @param[in]  level  The level
         *
         */
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


            this->get_multilevel_gradient(fun, x, g, g_diff, H_diff, s_global, level); 
            energy_old = this->get_multilevel_energy(fun,  x,  g_diff,  H_diff, s_global, level); 
            g_norm = norm2(g); 
            
            if(this->verbose() &&  verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE)
            {
                this->print_level_info(level); 
                PrintInfo::print_iter_status(0, {g_norm, energy_old, ared, pred, rho, get_delta(level-1) }); 
            }

            it++; 

            while(!converged)
            {
                
                this->get_multilevel_hessian(fun, x, H, H_diff, level); 
                energy_old = this->get_multilevel_energy(fun,  x,  g_diff,  H_diff, s_global, level); 

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

                energy_new = this->get_multilevel_energy(fun,  tp,  g_diff,  H_diff, s_global, level); 
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
                    this->get_multilevel_gradient(fun, x, g, g_diff, H_diff, s_global, level); 
                    g_norm = norm2(g); 
                    make_hess_updates =  this->update_hessian(g, g_old, s, H, rho, g_norm); 
                }

                converged  = (converged  == true) ? true : this->check_local_convergence(it_success,  g_norm, level, get_delta(level-1)); 
                
                if(this->verbose() &&  verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE)
                    PrintInfo::print_iter_status(it, {g_norm, energy_new, ared, pred, rho, get_delta(level-1)}); 
                it++; 

            }

            if(this->verbose() &&  verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE)
            {
                ColorModifier color_def(FG_DEFAULT);
                std::cout<< color_def; 
            }

            return reduction; 
        }





    protected:
        /**
         * @brief      Updates delta on given level 
         *
         * @param[in]  rho        The rho
         * @param[in]  level      The level
         * @param[in]  s_global   Sum of all corrections on given level 
         * @param      converged  convergence flag
         */
        virtual void delta_update(const Scalar & rho, const SizeType & level, const Vector & s_global, bool & converged)
        {
            Scalar intermediate_delta; 

            if(rho < this->eta1())
            {   
                 intermediate_delta = std::max(this->gamma1() * this->get_delta(level-1), 1e-15); 
            }
            else if (rho > this->eta2() )
            {
                 intermediate_delta = std::min(this->gamma2() * this->get_delta(level-1), 1e15); 
            }
            else
            {
                intermediate_delta = this->get_delta(level-1); 
            }      


            // on the finest level we work just with one radius 
            if(level==this->num_levels())
            {
                this->set_delta(level-1, intermediate_delta); 
            }
            else
            {
                Scalar corr_norm = this->level_dependent_norm(s_global, level); 
                converged = this->delta_termination(corr_norm, level); 
                if(converged)
                {
                    std::cout<<"termination  due to small radius for given level ... \n"; 
                    return; 
                }

                corr_norm = this->get_delta_zero(level-1) - corr_norm; 
                corr_norm = std::min(intermediate_delta, corr_norm); 

                if(corr_norm <= 0)
                    corr_norm = 0; 

                this->set_delta(level-1, std::min(intermediate_delta, corr_norm)); 
            }
        }



        /**
         * @brief      Checks for termination
         *
         * @param[in]  it_success  Number of succesfull iterations on given level
         * @param[in]  g_norm      Norm of gradient 
         * @param[in]  level       The level
         * @param[in]  delta       The delta
         *
         */
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

            return this->criticality_measure_termination(g_norm); 
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



        /**
         * @brief      Measures, how smooth is gradient 
         * @note       Needs to be overriden in case of $\infty$ norm, as norm of grad is not relevant anymore 
         *
         * @param[in]  g_norm  Norm of gradient 
         *
         */
        virtual bool criticality_measure_termination(const Scalar & g_norm)
        {
            // this should be fancier based on Graffon paper, but it is quite boring to work with so many mesh informations
            // Scalar _eps_grad_termination = std::min(0.001, eps_grad_termination[level]/ psi); 
            return (g_norm < _eps_grad_termination) ? true : false;    
        }


        /**
         * @brief      "Heuristics", which decides if it makes sense to go to the coarse level or no
         *
         * @param[in]  g_restricted  Restricted gradient
         * @param[in]  g_coarse      Coarse level gradient 
         *
         */
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


        /**
         * @brief      "Heuristics", which decides if hessian needs to be updated or now
         *
         * @param[in]  g_new   new gradient
         * @param[in]  g_old   Old gradient
         * @param[in]  s       The correction
         * @param[in]  H       The hessian 
         * @param[in]  rho     The rho
         * @param[in]  g_norm  Norm of gradient 
         *
         */
        virtual bool update_hessian(const Vector & g_new, const Vector & g_old, const Vector & s, const Matrix & H, const Scalar & rho, const Scalar & g_norm)
        {
            // iteration is not sucessful enough
            if(rho > 0 && rho < _hessian_update_eta)
                return true; 

            Vector help = g_new - g_old - H * s; 

            // Hessian approx is relativelly poor
            return (norm2(help) > _hessian_update_delta * g_norm) ? true : false; 
        }



        /**
         * @brief      Sets the delta.
         *
         * @param[in]  level   The level
         * @param[in]  radius  The radius
         *
         */
        virtual bool set_delta(const SizeType & level, const Scalar & radius)
        {
            _deltas[level] = radius; 
            return true; 
        }


        /**
         * @brief      Gets the delta.
         *
         * @param[in]  level  The level
         *
         * @return     The delta.
         */
        virtual Scalar get_delta(const SizeType & level) const 
        {
            return _deltas[level]; 
        }


        /**
         * @brief      Sets the delta zero.
         *
         * @param[in]  level   The level
         * @param[in]  radius  The radius
         *
         */
        virtual bool set_delta_zero(const SizeType & level, const Scalar & radius)
        {
            _deltas_zero[level] = radius; 
            return true; 
        }


        /**
         * @brief      Gets the delta zero.
         *
         * @param[in]  level  The level
         *
         * @return     The delta zero.
         */
        virtual Scalar get_delta_zero(const SizeType & level) const 
        {
            return _deltas_zero[level]; 
        }


        /**
         * @brief      Initializes tr radius on eaxh level. Organized from coarsest => delta[0] =  coarsest level
         *
         */
        virtual bool init_deltas()
        {
            for(Scalar i = 0; i < this->num_levels(); i ++)
                _deltas.push_back(_delta_init); 

            for(Scalar i = 0; i < this->num_levels()-1; i ++)
                _deltas_zero.push_back(0); 

            return true; 
        }


    
        /**
        * @brief      Initializes _delta_hessians for all levels. They are organized from coarsest to finest =>  _delta_hessians[0] =  coarsest level
         * @NOTE:      We do not have any _delta_hessians for the finest level, since function on the finest level is taken from problem discretization 
         *
         */
        bool init_delta_hessians()
        {
            _delta_hessians.resize(this->num_levels()-1); 
            return true; 
        }


        /**
         * @brief      Sets the delta hessian.
         *
         * @param[in]  level   The level
         * @param[in]  H_diff  The h difference
         *
         */
        virtual bool set_delta_hessian(const SizeType & level, const Matrix & H_diff)
        {
            _delta_hessians[level] = H_diff; 
            return true; 
        }


        /**
         * @brief      Gets the delta hessian.
         *
         * @param[in]  level  The level
         *
         * @return     The delta hessian.
         */
        virtual Matrix & get_delta_hessian(const SizeType & level) 
        {
            return _delta_hessians[level]; 
        }


        /**
         * @brief      Initializes x0 for all levels. They are organized from coarsest to finest =>  _x_initials[0] =  coarsest level.
         * @NOTE:      We do not have any x0 for the finest level, since function on the finest level is taken from problem discretization 
         *
         */
        virtual bool init_x_initials()
        {
            _x_initials.resize(this->num_levels()-1); 
            return true; 
        }


        /**
         * @brief      Sets the x initial.
         *
         * @param[in]  level  The level
         * @param[in]  x      The initial x.
         *
         */
        virtual bool set_x_initial(const SizeType & level, const Vector & x)
        {
            _x_initials[level] = x; 
            return true; 
        }


        /**
         * @brief      Gets the x initial.
         *
         * @param[in]  level  The level
         *
         * @return     The initial x.
         */
        virtual Vector & get_x_initial(const SizeType & level) 
        {
            return _x_initials[level]; 
        }

        
        /**
         * @brief      Initializes delta_grads for all levels. They are organized from coarsest to finest => _delta_gradients[0] =  coarsest level. 
         * @NOTE:      We do not have any g_diff for the finest level, since function on the finest level is taken from problem discretization 
         *
         */
        virtual bool init_delta_gradients()
        {
            _delta_gradients.resize(this->num_levels()-1); 
            return true; 
        }


        /**
         * @brief      Sets the delta gradient.
         *
         * @param[in]  level   The level
         * @param[in]  g_diff  The g_diff
         *
         */
        virtual bool set_delta_gradient(const SizeType & level, const Vector & g_diff)
        {
            _delta_gradients[level] = g_diff; 
            return true; 
        }


        /**
         * @brief      Gets the delta gradient.
         *
         * @param[in]  level  The level
         *
         * @return     The delta gradient.
         */
        virtual Vector & get_delta_gradient(const SizeType & level) 
        {
            return _delta_gradients[level]; 
        }


        

        /**
         * @brief      Provides coarse solve 
         *              it is here, just to be able to use full cycle 
         *              TODO: CHECK IF RHS does not need to be set-up
         *
         * @param      fun          The fun
         * @param      x            THe current iterate
         * @param[in]  rhs          The rhs
         *
         */
        virtual bool coarse_solve(FunctionType &fun, Vector &x, const Vector & /*rhs*/) override
        {
            local_tr_solve(fun, x, 0); 
            return true; 
        }


        /**
         * @brief      Prints some info related to level 
         *
         * @param[in]  level  The level
         */
        virtual void print_level_info(const SizeType & level)
        {
            ColorModifier color_out(FG_LIGHT_YELLOW);
            if(this->verbose() &&  verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE)
            {
                if(level == 1)
                {
                    std::cout << color_out; 
                    std::string solver_type = "COARSE SOLVE:: " + std::to_string(level); 
                    this->init_solver(solver_type, {" it. ", "|| g_norm ||", "   E + <g_diff, s>", "ared   ",  "  pred  ", "  rho  ", "  delta "}); 
                }
                else
                {
                    color_out.set_color_code(FG_LIGHT_GREEN); 
                    std::cout << color_out; 
                    std::string solver_type = "SMOOTHER:  " + std::to_string(level); 
                    this->init_solver(solver_type, {" it. ", "|| g_norm ||", "   E + <g_diff, s>", "ared   ",  "  pred  ", "  rho  ", "  delta "}); 
                }
            }
        }



        /**
         * @brief      Solves TR subroblem for given level 
         *
         * @param[in]  H      The hessian
         * @param[in]  g      The gradient
         * @param      s      New correction
         * @param[in]  level  The level
         *
         */
        virtual bool solve_qp_subproblem(const Matrix & H, const Vector & g, Vector & s, const SizeType & level)
        {
            if(level == 1)
            {
                _coarse_tr_subproblem->current_radius(get_delta(level-1));  
                _coarse_tr_subproblem->atol(1e-16); 
                _coarse_tr_subproblem->tr_constrained_solve(H, g, s); 
            }
            else
            {
                _smoother_tr_subproblem->current_radius(get_delta(level-1));  
                _smoother_tr_subproblem->atol(1e-16); 
                _smoother_tr_subproblem->max_it(5);
                _smoother_tr_subproblem->tr_constrained_solve(H, g, s); 

            }

            return true; 
        }


        /**
         * @brief      Computes hessian on given level 
         *
         * @param[in]  fun       Function with evaluation routines
         * @param[in]  x         The current iterate
         * @param      H         The current hessian
         * @param[in]  H_diff    The h_diff
         * @param[in]  level     The level
         *
         * @return     The multilevel hessian.
         */
        virtual bool get_multilevel_hessian(const FunctionType & fun, const Vector & x,  Matrix & H, const Matrix & H_diff, const SizeType & level)
        {
            if(level < this->num_levels())
                return MultilevelHessianEval<Matrix, Vector, FunctionType, CONSISTENCY_LEVEL>::compute_hessian(fun, x, H, H_diff);
            else
                return fun.hessian(x, H); 
        }


        /**
         * @brief      Computes gradient on given level 
         *
         * @param[in]  fun       Function with evaluation routines
         * @param[in]  x         The current iterate
         * @param      g         The current gradient
         * @param[in]  g_diff    The g_diff
         * @param[in]  H_diff    The h_diff
         * @param[in]  s_global  The sum of all corrections on given level 
         * @param[in]  level     The level
         *
         * @return     The multilevel gradient.
         */
        virtual bool get_multilevel_gradient(const FunctionType & fun, const Vector & x,  Vector & g, const Vector & g_diff, const Matrix & H_diff, const Vector & s_global, const SizeType & level)
        {
            if(level < this->num_levels())
                return MultilevelGradientEval<Matrix, Vector, FunctionType, CONSISTENCY_LEVEL>::compute_gradient(fun, x, g, g_diff, H_diff, s_global);
            else
                 return fun.gradient(x, g); 
        }

        /**
         * @brief      Computes energy for given level 
         *
         * @param[in]  fun       Function with evaluation routines
         * @param[in]  x         The current iterate
         * @param[in]  g_diff    The g_diff
         * @param[in]  H_diff    The h_diff
         * @param[in]  s_global  The sum of all corrections on given level 
         * @param[in]  level     The level
         *
         * @return     The multilevel energy.
         */
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


        /**
         * @brief      Computes norm of coarse level vector wrt to fine level
         *
         * @param[in]  u          The current iterate
         * @param[in]  current_l  The current level
         *
         */
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




        virtual void print_statistics() const 
        {
            auto rmtr_data_path = Utopia::Instance().get("rmtr_data_path");
            if(!rmtr_data_path.empty())
            {
                CSVWriter writer; 
                if(!writer.file_exists(rmtr_data_path))
                {
                    writer.open_file(rmtr_data_path); 
                    writer.write_table_row<std::string>({("v_cycles")}); 
                }
                else
                    writer.open_file(rmtr_data_path); 

                writer.write_table_row<SizeType>({(_it_global)}); 
                writer.close_file(); 
            }
        }




    protected:   
        SizeType                            _it_global;                 /** * global iterate counter  */
        std::vector<Scalar>                 _deltas;                    /** * deltas on given level  */
        std::vector<Scalar>                 _deltas_zero;               /** * initial deltas on given level  */

        


        std::shared_ptr<TRSubproblem>        _coarse_tr_subproblem;     /** * solver used to solve coarse level TR subproblems  */
        std::shared_ptr<TRSubproblem>        _smoother_tr_subproblem;   /** * solver used to solve fine level TR subproblems  */


        std::vector<Vector>           _delta_gradients;             /** * difference between fine and coarse level gradient */
        std::vector<Matrix>           _delta_hessians;              /** * difference between fine and coarse level hessians */
        std::vector<Vector>           _x_initials;                  /** * initial iterates on given level */


        // ----------------------- PARAMETERS ----------------------
        Parameters                      _parameters; 

        
        Scalar                          _delta_init;                /** * delta zero  */
        SizeType                        _max_coarse_it;             /** * maximum iterations on coarse level   */
        SizeType                        _max_smoothing_it;          /** * max smoothing iterations  */

        Scalar                         _eps_delta_termination;      /** * maximum delta allowed on coarse level - makes sure that coarse level corection stays inside fine level radius  */
        Scalar                         _delta_min;                  /** * minimum delta allowed by algorithm   */

        Scalar                         _grad_smoothess_termination; /** * determines when gradient is not smooth enough => does pay off to go to coarse level at all  */
        Scalar                         _eps_grad_termination;       /** * tolerance on grad  */

        Scalar                         _hessian_update_delta;       /** * tolerance used for updating hessians */
        Scalar                         _hessian_update_eta;         /** * tolerance used for updating hessians */


        VerbosityLevel                  _verbosity_level; 


    };

}

#endif //UTOPIA_RMTR_HPP

