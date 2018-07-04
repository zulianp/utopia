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


#include "utopia_LevelMemory.hpp"

namespace utopia 
{
    /**
     * @brief      The class for Recursive multilevel trust region solver. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
    class RMTR : public NonlinearMultiLevelBase<Matrix, Vector>,
                       public TrustRegionBase<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;
        typedef utopia::TRSubproblem<Matrix, Vector>        TRSubproblem; 
        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef utopia::Level<Matrix, Vector>               Level;
        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;

    public:

        using TrustRegionBase<Matrix, Vector>::delta_update;
        using NonlinearMultiLevelBase<Matrix, Vector>::print_statistics;

       /**
        * @brief      Multigrid class
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level. 
        */
        RMTR(   const std::shared_ptr<TRSubproblem> &tr_subproblem_coarse,  const std::shared_ptr<TRSubproblem> &tr_subproblem_smoother, 
                const Parameters params = Parameters()): 
                NonlinearMultiLevelBase<Matrix,Vector>(params), 
                _coarse_tr_subproblem(tr_subproblem_coarse), 
                _smoother_tr_subproblem(tr_subproblem_smoother) 
        {
            set_parameters(params); 
        }

        virtual ~RMTR(){} 
        

        void set_parameters(const Parameters params) override
        {
            NonlinearMultiLevelBase<Matrix, Vector>::set_parameters(params);    
            _it_global                  = 0;          
            _parameters                 = params; 

            _max_coarse_it              = params.max_coarse_it();
            _max_smoothing_it           = params.max_smoothing_it();
            _eps_delta_termination      = params.eps_delta_termination();
            _grad_smoothess_termination = params.grad_smoothess_termination();
            _eps_grad_termination       = params.eps_grad_termination();
            _hessian_update_delta       = params.hessian_update_delta();
            _hessian_update_eta         = params.hessian_update_eta();

            _verbosity_level           = params.verbosity_level(); 
        }

        virtual void init_memory(const SizeType & fine_local_size) override 
        {
            std::cout<<"-------- to be done \n"; 

            memory_.init(this->n_levels()); 
            memory_.g_diff[this->n_levels()-1] = local_zeros(fine_local_size); 
        }



        VerbosityLevel verbosity_level() const 
        {
            return _verbosity_level; 
        }

        void verbosity_level(const VerbosityLevel & level )
        {

            _verbosity_level = this->verbose() ? level : VERBOSITY_LEVEL_QUIET;  
        }


        void set_grad_smoothess_termination(const Scalar & grad_smoothess_termination)
        {
            _grad_smoothess_termination = grad_smoothess_termination; 
        }


        using NonlinearMultiLevelBase<Matrix, Vector>::solve; 

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
        virtual bool solve(Fun &fine_fun, Vector & x_h, const Vector & rhs) override
        {
            bool converged = false; 
            SizeType l = this->n_levels(); 
            Scalar r_norm, r0_norm, rel_norm, energy;

            Vector g_finest  = local_zeros(local_size(x_h)); 
            this->make_iterate_feasible(fine_fun, x_h); 

            fine_fun.gradient(x_h, g_finest); 
            fine_fun.value(x_h, energy); 

            r0_norm = norm2(g_finest); 
            _it_global = 0; 

            //-------------- INITIALIZATIONS ---------------
            this->status_.clear();
            init(); 
            
            //----------------------------------------------

            if(verbosity_level() >= VERBOSITY_LEVEL_NORMAL && mpi_world_rank() == 0)
            {
                ColorModifier red(FG_LIGHT_MAGENTA);
                ColorModifier def(FG_DEFAULT);
                std::cout << red;

                std::string name_id = this->name_id() + "     Number of levels: " + std::to_string(l); 
                this->init_solver(name_id, {" it. ", "|| g_norm ||", "   E "}); 

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

                fine_fun.gradient(x_h, g_finest); 
                fine_fun.value(x_h, energy); 
                
                r_norm = norm2(g_finest);
                rel_norm = r_norm/r0_norm; 

                _it_global++; 

                if(this->verbose() && mpi_world_rank() == 0)
                {
                    ColorModifier red(FG_LIGHT_MAGENTA);
                    ColorModifier def(FG_DEFAULT);
                    std::cout << red; 

                    if(verbosity_level() > VERBOSITY_LEVEL_NORMAL)
                        this->print_init_message("RMTR OUTER SOLVE", {" it. ", "|| g_norm ||", "   E "}); 

                    PrintInfo::print_iter_status(_it_global, {r_norm, energy}); 
                    std::cout << def; 
                }

                // check convergence
                converged = check_global_convergence(_it_global, r_norm, rel_norm, this->get_delta(l-1)); 
            }

            // benchmarking
            print_statistics(); 
            return true; 
        }



    private: 

        /**
         * @brief      Multiplicative cycle 
         *
         * @param      fine_fun   The fine fun
         * @param      u_l        Current iterate
         * @param[in]  f          The right hand side
         * @param[in]  level      The level
         *
         */
        virtual bool multiplicative_cycle(Fun &fine_fun, Vector & u_l, const Vector &/*f*/, const SizeType & level) override
        {
            Vector g_fine, g_coarse, g_diff, u_2l, s_coarse, s_fine, s_global; 
            Matrix H_fine, H_coarse, H_diff; 

            Scalar ared=0.0, coarse_reduction=0.0, rho=0.0; 
            Scalar E_old, E_new; 
            bool converged = false, smoothness_flg=true; 

            //----------------------------------------------------------------------------
            //                   presmoothing
            //----------------------------------------------------------------------------
            converged = this->local_tr_solve(fine_fun, u_l, level); 

            // making sure that correction does not exceed tr radius ... 
            if(converged)
                return true; 

            compute_s_global(u_l, level, s_global); 
            this->get_multilevel_gradient(fine_fun, u_l, g_fine, s_global, level); 

            if(level == this->n_levels())
            {
                converged =  this->criticality_measure_termination(norm2(g_fine)); 
                if(converged==true)
                    return true; 
            }

            this->transfer(level-2).restrict(g_fine, g_diff);
            this->transfer(level-2).project_down(u_l, u_2l); 

            this->make_iterate_feasible(this->function(level-2), u_2l); 

            //----------------------------------------------------------------------------
            //                   first order coarse level objective managment
            //----------------------------------------------------------------------------            
            if(CONSISTENCY_LEVEL != GALERKIN)
            {             
                this->function(level-2).gradient(u_2l, g_coarse); 
                this->zero_correction_related_to_equality_constrain(this->function(level-2), g_diff); 
            }

            smoothness_flg = grad_smoothess_termination(g_diff, g_fine); 

            if(CONSISTENCY_LEVEL != GALERKIN)
                g_diff -= g_coarse; 

            //----------------------------------------------------------------------------
            //                   second order coarse level objective managment
            //----------------------------------------------------------------------------            
            if(CONSISTENCY_LEVEL == SECOND_ORDER || CONSISTENCY_LEVEL == GALERKIN)
            {
                this->get_multilevel_hessian(fine_fun, u_l, H_fine, level); 
                this->transfer(level-2).restrict(H_fine, H_diff);
                
                if(CONSISTENCY_LEVEL == SECOND_ORDER)
                {
                    this->zero_correction_related_to_equality_constrain_mat(this->function(level-2), H_diff); 
                    this->function(level-2).hessian(u_2l, H_coarse); 
                    H_diff -=  H_coarse; 
                }
            }


            //----------------------------------------------------------------------------
            //                   initializing coarse level
            //----------------------------------------------------------------------------
            this->set_delta(level-2, get_delta(level-1)); 

            this->set_delta_gradient(level-2, g_diff); 
            if(CONSISTENCY_LEVEL == SECOND_ORDER || CONSISTENCY_LEVEL == GALERKIN)
                this->set_delta_hessian(level-2, H_diff); 

            this->set_x_initial(level-2, u_2l); 
            
        
            s_coarse = 0*u_2l; 
            coarse_reduction = this->get_multilevel_energy(this->function(level-2),  u_2l,  s_coarse, level-1); 


            //----------------------------------------------------------------------------
            //               recursion  / Taylor correction
            //----------------------------------------------------------------------------
            if(level == 2 && smoothness_flg)
            {
                SizeType l_new = level - 1; 
                this->local_tr_solve(this->function(level-2), u_2l, l_new, true); 
            }
            else if(smoothness_flg)
            {
                // recursive call into RMTR 
                for(SizeType k = 0; k < this->mg_type(); k++)
                {   
                    SizeType l_new = level - 1; 
                    this->multiplicative_cycle(this->function(level-2), u_2l, g_diff, l_new); 
                }
            }
            
            if(smoothness_flg)
            {                
                //----------------------------------------------------------------------------
                //                       building trial point from coarse level 
                //----------------------------------------------------------------------------
                s_coarse = u_2l - this->get_x_initial(level - 2);
                coarse_reduction -= this->get_multilevel_energy(this->function(level-2),  u_2l,  s_coarse, level-1);

                this->transfer(level-2).interpolate(s_coarse, s_fine);
                this->zero_correction_related_to_equality_constrain(fine_fun, s_fine); 

                compute_s_global(u_l, level, s_global);                               
                E_old = this->get_multilevel_energy(fine_fun,  u_l,  s_global, level); 

                // new test for dbg mode 
                u_l += s_fine; 

                compute_s_global(u_l, level, s_global);     
                E_new = this->get_multilevel_energy(fine_fun,  u_l,  s_global, level); 
                
                //----------------------------------------------------------------------------
                //                        trial point acceptance  
                //----------------------------------------------------------------------------

                ared = E_old - E_new; 
                rho = ared / coarse_reduction; 
                if(coarse_reduction<=0)
                    rho = 0; 

                Scalar coarse_corr_taken = 0; 
                if(rho > this->rho_tol())
                {
                    coarse_corr_taken = 1; 
                }
                else
                {
                    u_l -= s_fine; 
                    compute_s_global(u_l, level, s_global);
                }

                //----------------------------------------------------------------------------
                //                                  trust region update 
                //----------------------------------------------------------------------------
                converged = delta_update(rho, level, s_global); 
                
                // terminate, since TR rad. does not allow to take more corrections on given level 
                if(converged==true) 
                    return true; 

                if(verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
                {
                    // just to see what is being printed 
                    std::string status = "RMTR_coarse_corr_stat, level: " + std::to_string(level); 
                    this->print_init_message(status, {" it. ", "   E_old     ", "   E_new", "ared   ",  "  coarse_level_reduction  ", "  rho  ", "  delta ", "taken"}); 
                    PrintInfo::print_iter_status(_it_global, {E_old, E_new, ared, coarse_reduction, rho, get_delta(level-1), coarse_corr_taken }); 
                }
            }
            else
            {
                std::cout<< "-------------------------------------- GRAD NOT SMOOTH ENOUGH ----------------------------------- \n"; 
            }

            //----------------------------------------------------------------------------
            //                        postsmoothing   
            //----------------------------------------------------------------------------
            this->local_tr_solve(fine_fun, u_l, level, !smoothness_flg); 
    
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
        virtual bool local_tr_solve(Fun &fun, Vector & x, const SizeType & level, const bool & exact_solve_flg = false)
        {   
            Vector s_global, g; 
            Matrix  H; 

            SizeType it_success = 0, it = 0; 
            Scalar ared = 0. , pred = 0., rho = 0., energy_old=9e9, energy_new=9e9, g_norm=1.0; 
            bool make_grad_updates = true, /*make_hess_updates = true,*/ converged = false, delta_converged = false; 

            Vector s = local_zeros(local_size(x)); 
            
            compute_s_global(x, level, s_global);  
            this->get_multilevel_gradient(fun, x, g, s_global, level); 
            energy_old = this->get_multilevel_energy(fun,  x, s_global, level); 
            g_norm = norm2(g); 


            if(verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
            {
                this->print_level_info(level); 
                PrintInfo::print_iter_status(0, {g_norm, energy_old, ared, pred, rho, get_delta(level-1) }); 
            }

            it++;       

            while(!converged)
            {
                this->get_multilevel_hessian(fun, x, H, level); 

            //----------------------------------------------------------------------------
            //     solving constrained system to get correction and  building trial point 
            //----------------------------------------------------------------------------
                // correction needs to get prepared 
                s = 0 * x;
                this->solve_qp_subproblem(H, g, s, level, exact_solve_flg); 

                // predicted reduction based on model 
                TrustRegionBase<Matrix, Vector>::get_pred(g, H, s, pred); 

                // building trial point 
                x += s;  
            

                compute_s_global(x, level, s_global); 
                energy_new = this->get_multilevel_energy(fun,  x, s_global, level); 
                ared = energy_old - energy_new; 
                
                rho = (ared < 0) ? 0.0 : ared/pred; 
                rho = (rho != rho) ? 0.0 : rho; 

            //----------------------------------------------------------------------------
            //     acceptance of trial point 
            //----------------------------------------------------------------------------
              
                // good reduction, accept trial point 
                if (rho >= this->rho_tol())
                {
                    it_success++; 
                    make_grad_updates =  true; 
                    energy_old = energy_new; 
                }
                else
                {   
                    x -= s; // return iterate into initial state 
                    compute_s_global(x, level, s_global); 
                    make_grad_updates =  false; 
                }
            //----------------------------------------------------------------------------
            //     trust region update 
            //----------------------------------------------------------------------------
                delta_converged = delta_update(rho, level, s_global); 

                if(make_grad_updates)
                {
                    Vector g_old = g; 
                    this->get_multilevel_gradient(fun, x, g, s_global, level); 
                    g_norm = norm2(g); 
                    // make_hess_updates =  
                    this->update_hessian(g, g_old, s, H, rho, g_norm); 
                }

                if(verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
                    PrintInfo::print_iter_status(it, {g_norm, energy_new, ared, pred, rho, get_delta(level-1)}); 

                converged  = (delta_converged  == true) ? true : this->check_local_convergence(it_success,  g_norm, level, get_delta(level-1)); 
                
                if(level == this->n_levels())
                    converged  = (converged  == true || g_norm < this->atol_) ? true : false; 

                it++; 

            }

            if(verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
            {
                ColorModifier color_def(FG_DEFAULT);
                std::cout<< color_def; 
            }

            return delta_converged; 
        }



    protected:


        virtual void init() 
        {
            // new version .....
            memory_.init(this->n_levels()); 
            
            // new init deltas... 
            for(Scalar l = 0; l < this->n_levels(); l ++)
                memory_.delta[l] = this->delta0(); 
        }

        // -------------------------- tr radius managment ---------------------------------------------        
        /**
         * @brief      Updates delta on given level 
         *
         * @param[in]  rho        The rho
         * @param[in]  level      The level
         * @param[in]  s_global   Sum of all corrections on given level 
         * @param      converged  convergence flag
         */
        virtual bool delta_update(const Scalar & rho, const SizeType & level, const Vector & s_global)
        {
            Scalar intermediate_delta; 

            if(rho < this->eta1())
                 intermediate_delta = std::max(this->gamma1() * this->get_delta(level-1), 1e-15); 
            else if (rho > this->eta2() )
                 intermediate_delta = std::min(this->gamma2() * this->get_delta(level-1), 1e15); 
            else
                intermediate_delta = this->get_delta(level-1); 


            // on the finest level we work just with one radius 
            if(level==this->n_levels())
            {
                this->set_delta(level-1, intermediate_delta); 
                return false; 
            }
            else
            {
                Scalar corr_norm = this->level_dependent_norm(s_global, level); 
                bool converged = this->delta_termination(corr_norm, level); 
                
                if(converged && verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
                    std::cout<<"termination  due to small radius on level: "<< level << ". \n"; 

                corr_norm = this->get_delta(level) - corr_norm; 
                corr_norm = std::min(intermediate_delta, corr_norm); 

                if(corr_norm <= 0.0)
                    corr_norm = 0.0; 

                this->set_delta(level-1, corr_norm); 
                return converged; 
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
            if(current_l == this->n_levels())
                return 0.0; 
            else
            {
                Vector s = u; // carries over prolongated correction
                for(SizeType i = current_l; i < this->n_levels(); i++)
                    this->transfer(i-1).interpolate(s, s); 
                
                return norm2(s); 
            }    
        }



    // ---------------------------------- convergence checks -------------------------------
        virtual bool check_global_convergence(const SizeType & it, const Scalar & r_norm, const Scalar & rel_norm, const Scalar & delta)
        {   
            bool converged = NonlinearMultiLevelBase<Matrix, Vector>::check_convergence(it, r_norm, rel_norm, 1); 

            if(delta < this->delta_min())
            {
                converged = true; 
                this->exit_solver(it, ConvergenceReason::CONVERGED_TR_DELTA);
            }

            return converged; 
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

            if(delta < this->delta_min())
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
            return (corr_norm > (1.0 - _eps_delta_termination) * get_delta(level)) ? true : false; 
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
            return (Rg_norm >= _grad_smoothess_termination * g_norm) ? true : false;   
        }


// ------------------------------- efficiency  thing --------------------
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



// ------------------------------- initializations  --------------------

        /**
         * @brief      Sets the delta.
         *
         * @param[in]  level   The level
         * @param[in]  radius  The radius
         *
         */
        virtual bool set_delta(const SizeType & level, const Scalar & radius)
        {
            memory_.delta[level] = radius; 
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
            return memory_.delta[level]; 
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
            memory_.H_diff[level] = H_diff; 
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
            return memory_.H_diff[level]; 
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
            memory_.x_0[level] = x; 
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
            return memory_.x_0[level];
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
            memory_.g_diff[level] = g_diff; 
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
            return memory_.g_diff[level]; 
        }


//----------------------------- QP solve -----------------------------------------------------------------


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
        virtual bool coarse_solve(Fun &fun, Vector &x, const Vector & /*rhs*/) override
        {
            local_tr_solve(fun, x, 0); 
            return true; 
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
        virtual bool solve_qp_subproblem(const Matrix & H, const Vector & g, Vector & s, const SizeType & level, const bool & flg)
        {
            if(flg)
            {
                _coarse_tr_subproblem->atol(1e-16); 
                _coarse_tr_subproblem->max_it(5000); 
                _coarse_tr_subproblem->tr_constrained_solve(H, g, s, get_delta(level-1)); 
            }
            else
            {
                _smoother_tr_subproblem->atol(1e-16); 
                _smoother_tr_subproblem->max_it(5);
                _smoother_tr_subproblem->tr_constrained_solve(H, g, s, get_delta(level-1)); 
            }

            return true; 
        }


//----------------------------------------- energy evaluation helpers -------------------------------------------

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
        virtual bool get_multilevel_hessian(const Fun & fun, const Vector & x,  Matrix & H, const SizeType & level)
        {
            if(level < this->n_levels())
                return MultilevelHessianEval<Matrix, Vector, CONSISTENCY_LEVEL>::compute_hessian(fun, x, H, get_delta_hessian(level-1));
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
        virtual bool get_multilevel_gradient(const Fun & fun, const Vector & x,  Vector & g, const Vector & s_global, const SizeType & level)
        {
            if(level < this->n_levels())
            {
                return MultilevelGradientEval<Matrix, Vector, CONSISTENCY_LEVEL>::compute_gradient(fun, x, g, get_delta_gradient(level-1), get_delta_hessian(level-1), s_global);
            }
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
        virtual Scalar get_multilevel_energy(const Fun & fun, const Vector & x, const Vector & s_global, const SizeType & level) 
        {
            if(level < this->n_levels())
                return MultilevelEnergyEval<Matrix, Vector, CONSISTENCY_LEVEL>::compute_energy(fun, x, get_delta_gradient(level-1), get_delta_hessian(level-1), s_global); 
            else
            {
                Scalar energy; 
                fun.value(x, energy); 
                return energy; 
            }
        }


        virtual void compute_s_global(const Vector & x, const SizeType & level, Vector & s_global)
        {
            if(level < this->n_levels())
                s_global = x - this->get_x_initial(level - 1);         
        }


//----------------------------------------- printouts helpers -------------------------------------------

        /**
         * @brief      Prints some info related to level 
         *
         * @param[in]  level  The level
         */
        virtual void print_level_info(const SizeType & level)
        {
            ColorModifier color_out(FG_LIGHT_YELLOW);
            if(verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
            {
                if(level == 1)
                {
                    std::cout << color_out; 
                    std::string solver_type = "COARSE SOLVE:: " + std::to_string(level); 
                    this->print_init_message(solver_type, {" it. ", "|| g_norm ||", "   E + <g_diff, s>", "ared   ",  "  pred  ", "  rho  ", "  delta "}); 
                }
                else
                {
                    color_out.set_color_code(FG_LIGHT_GREEN); 
                    std::cout << color_out; 
                    std::string solver_type = "SMOOTHER:  " + std::to_string(level); 
                    this->print_init_message(solver_type, {" it. ", "|| g_norm ||", "   E + <g_diff, s>", "ared   ",  "  pred  ", "  rho  ", "  delta "}); 
                }
            }
        }



        virtual void print_statistics() 
        {
            auto rmtr_data_path = Utopia::instance().get("rmtr_data_path");
            if(!rmtr_data_path.empty())
            {
                CSVWriter writer; 
                if (mpi_world_rank() == 0)
                {
                    if(!writer.file_exists(rmtr_data_path))
                    {
                        writer.open_file(rmtr_data_path); 
                        writer.write_table_row<std::string>({"v_cycles", "time"}); 
                    }
                    else
                        writer.open_file(rmtr_data_path); 

                    writer.write_table_row<Scalar>({Scalar(_it_global), this->get_time()}); 
                    writer.close_file(); 
                }
            }
        }




    protected:   
        SizeType                            _it_global;                 /** * global iterate counter  */

        std::shared_ptr<TRSubproblem>        _coarse_tr_subproblem;     /** * solver used to solve coarse level TR subproblems  */
        std::shared_ptr<TRSubproblem>        _smoother_tr_subproblem;   /** * solver used to solve fine level TR subproblems  */


        // ----------------------- PARAMETERS ----------------------
        Parameters                      _parameters; 

    
        SizeType                        _max_coarse_it;             /** * maximum iterations on coarse level   */
        SizeType                        _max_smoothing_it;          /** * max smoothing iterations  */

        Scalar                         _eps_delta_termination;      /** * maximum delta allowed on coarse level - makes sure that coarse level corection stays inside fine level radius  */

        Scalar                         _grad_smoothess_termination; /** * determines when gradient is not smooth enough => does pay off to go to coarse level at all  */
        Scalar                         _eps_grad_termination;       /** * tolerance on grad  */

        Scalar                         _hessian_update_delta;       /** * tolerance used for updating hessians */
        Scalar                         _hessian_update_eta;         /** * tolerance used for updating hessians */


        VerbosityLevel                  _verbosity_level; 



    private:
        RMTRLevelMemory <Matrix, Vector>         memory_;


    };

}

#endif //UTOPIA_RMTR_HPP