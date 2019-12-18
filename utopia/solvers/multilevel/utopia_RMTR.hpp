#ifndef UTOPIA_RMTR_HPP
#define UTOPIA_RMTR_HPP
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"
#include "utopia_RMTRBase.hpp"

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

    enum LocalSolveType{    PRE_SMOOTHING  = 1,
                            POST_SMOOTHING = 2,
                            COARSE_SOLVE   = 0};


    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
    class RMTR : public RMTRBase<Matrix, Vector, CONSISTENCY_LEVEL>,
                 public TrustRegionBase<Vector>
    {
        typedef UTOPIA_SCALAR(Vector)           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)        SizeType;

        typedef utopia::TRSubproblem<Matrix, Vector>    TRSubproblem;
        typedef std::shared_ptr<TRSubproblem>           TRSubproblemPtr;

        typedef utopia::Transfer<Matrix, Vector>                        Transfer;
        typedef utopia::Level<Matrix, Vector>                           Level;
        
        typedef utopia::RMTRBase<Matrix, Vector, CONSISTENCY_LEVEL>     RMTRBase;
        typedef typename RMTRBase::Fun Fun;

    public:
        using TrustRegionBase<Vector>::delta_update;
        using TrustRegionBase<Vector>::get_pred;

       /**
        * @brief      Implementation of RMTR with tr radius defined by l2 norm 
        *
        * @param[in]  n_levels       Number of levels
        */
        RMTR(   const SizeType & n_levels): RMTRBase(n_levels)
        {

        }

        virtual ~RMTR(){}

        virtual void read(Input &in) override
        {
            RMTRBase::read(in);
            TrustRegionBase<Vector>::read(in);
            const auto n_subproblems = _tr_subproblems.size();

            if(n_subproblems > 0)
            {
                in.get("coarse-QPSolver", *_tr_subproblems[0]);

                for(std::size_t i = 1; i < n_subproblems; i++)
                    in.get("fine-QPSolver", *_tr_subproblems[i]);
            }
        }

        virtual void print_usage(std::ostream &os) const override
        {
            RMTRBase::print_usage(os);
            TrustRegionBase<Vector>::print_usage(os);

            this->print_param_usage(os, "coarse-QPSolver", "TRSubproblem", "Input parameters for fine level QP solvers.", "-");
            this->print_param_usage(os, "fine-QPSolver", "TRSubproblem", "Input parameters for coarse level QP solver.", "-");
        }

        using NonlinearMultiLevelBase<Matrix, Vector>::solve;

        virtual std::string name() override { return "RMTR";  }


        bool set_coarse_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy)
        {
            if(static_cast<SizeType>(_tr_subproblems.size()) != this->n_levels()){
                _tr_subproblems.resize(this->n_levels());
            }

            _tr_subproblems[0] = strategy;

            return true;
        }

        bool set_fine_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy)
        {
            if(static_cast<SizeType>(_tr_subproblems.size()) != this->n_levels()){
                _tr_subproblems.resize(this->n_levels());
            }

            // starting from level 1 ....
            for(auto l = 1; l != this->n_levels(); ++l){
                _tr_subproblems[l] = std::shared_ptr<TRSubproblem>(strategy->clone());
            }

            return true;
        }


        bool set_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy, const SizeType & level)
        {
            if(static_cast<SizeType>(_tr_subproblems.size()) != this->n_levels()){
                _tr_subproblems.resize(this->n_levels());
            }

            if(level <= this->n_levels()){
                _tr_subproblems[level] = strategy;
            }
            else{
                utopia_error("utopia::RMTR::set_tr_strategy:: Requested level exceeds number of levels in ML hierarchy. \n");
            }

            return true;
        }


        bool set_tr_strategies(const std::vector<TRSubproblemPtr> &strategies)
        {
            if(static_cast<SizeType>(strategies.size()) != this->n_levels()){
                utopia_error("utopia::RMTR::set_tr_strategies:: Number of tr strategies MUST be equal to number of levels in ML hierarchy. \n");
            }

            _tr_subproblems = strategies;
            return true;
        }


    protected:

        /**
         * @brief      Multiplicative cycle
         *
         * @param      fine_fun   The fine fun
         * @param      u_l        Current iterate
         * @param[in]  f          The right hand side
         * @param[in]  level      The level
         *
         */
        virtual bool multiplicative_cycle(const SizeType & level) override
        {
            Scalar ared=0.0, coarse_reduction=0.0, rho=0.0;
            Scalar E_old=0.0, E_new=0.0;
            bool converged = false, smoothness_flg=true;

            //----------------------------------------------------------------------------
            //                   presmoothing
            //----------------------------------------------------------------------------
            if(this->pre_smoothing_steps()!=0){
                converged = this->local_tr_solve(level, PRE_SMOOTHING);
            }
            else
            {
                converged = false;
            }


            // making sure that correction does not exceed tr radius ...
            if(converged){
                return true;
            }

            if(this->pre_smoothing_steps()==0 && level < this->n_levels()-1)
            {
                this->compute_s_global(level, this->memory_.s_working[level]);
                this->get_multilevel_gradient(this->function(level), this->memory_.s_working[level], level);
            }

            if(level == this->n_levels()-1)
            {
                converged =  this->criticality_measure_termination(this->criticality_measure(level));
                if(converged==true)
                    return true;
            }

            // Restricted fine level gradient 
            this->transfer(level-1).restrict(this->ml_derivs_.g[level], this->ml_derivs_.g_diff[level-1]);

            // Projecting current iterate to obtain initial iterate on coarser grid 
            this->transfer(level-1).project_down(this->memory_.x[level], this->memory_.x[level-1]);


            if(!this->skip_BC_checks()){
                if(CONSISTENCY_LEVEL != GALERKIN){
                    this->make_iterate_feasible(this->function(level-1), this->memory_.x[level-1]);
                }
            }

            //----------------------------------------------------------------------------
            //    initializing coarse level (deltas, constraints, hessian approx, ...)
            //----------------------------------------------------------------------------
            this->init_level(level-1);

            //----------------------------------------------------------------------------
            //                  first order coarse level objective managment
            //----------------------------------------------------------------------------
            if(CONSISTENCY_LEVEL != GALERKIN)
            {
                if(empty(this->ml_derivs_.g[level-1])){
                    this->ml_derivs_.g[level-1] = 0.0* this->memory_.x[level-1]; 
                }

                this->function(level-1).gradient(this->memory_.x[level-1], this->ml_derivs_.g[level-1]);
            }

            if(!this->skip_BC_checks())
            {
                if(CONSISTENCY_LEVEL != GALERKIN){
                    this->zero_correction_related_to_equality_constrain(this->function(level-1), this->ml_derivs_.g_diff[level-1]);
                }
            }

            if(this->check_grad_smoothness() && CONSISTENCY_LEVEL != GALERKIN){
                smoothness_flg = this->grad_smoothess_termination(this->ml_derivs_.g_diff[level-1], this->ml_derivs_.g[level-1], level-1);
            }
            else{
                smoothness_flg = true;
            }


            if(CONSISTENCY_LEVEL != GALERKIN)
            {
                this->ml_derivs_.g_diff[level-1] -= this->ml_derivs_.g[level-1];
            }

            //----------------------------------------------------------------------------
            //                   second order coarse level objective managment
            //----------------------------------------------------------------------------
            if(CONSISTENCY_LEVEL == SECOND_ORDER || CONSISTENCY_LEVEL == GALERKIN)
            {
                this->get_multilevel_hessian(this->function(level), level);
                this->transfer(level-1).restrict(this->ml_derivs_.H[level], this->ml_derivs_.H_diff[level-1]);

                if(CONSISTENCY_LEVEL == SECOND_ORDER)
                {
                    if(!this->skip_BC_checks()){
                        this->zero_correction_related_to_equality_constrain_mat(this->function(level-1), this->ml_derivs_.H_diff[level-1]);
                    }

                    this->function(level-1).hessian(this->memory_.x[level-1], this->ml_derivs_.H[level-1]);

                    // memory_.H_diff[level-1] = memory_.H_diff[level-1] -  memory_.H[level-1];
                    this->ml_derivs_.H_diff[level-1] -= this->ml_derivs_.H[level-1];
                }
            }

            //----------------------------------------------------------------------------
            //                   additional coarse level initialization...
            //----------------------------------------------------------------------------
            this->memory_.x_0[level-1]    = this->memory_.x[level-1];
            // this->memory_.s[level-1]      = local_zeros(local_size(this->memory_.x[level-1]));
            // todo:: check if necessary
            this->memory_.s[level-1].set(0.0); 

            // at this point s_global on coarse level is empty
            coarse_reduction = this->get_multilevel_energy(this->function(level-1), this->memory_.s[level-1], level-1);

            //----------------------------------------------------------------------------
            //               recursion  / Taylor correction
            //----------------------------------------------------------------------------
            if(level == 1 && smoothness_flg)
            {
                this->local_tr_solve(level - 1, COARSE_SOLVE);
            }
            else if(smoothness_flg)
            {
                // recursive call into RMTR
                for(SizeType k = 0; k < this->mg_type(); k++)
                {
                    SizeType l_new = level - 1;
                    this->multiplicative_cycle(l_new);
                }
            }

            if(smoothness_flg)
            {
                //----------------------------------------------------------------------------
                //                       building trial point from coarse level
                //----------------------------------------------------------------------------
                this->memory_.s[level-1] = this->memory_.x[level-1] - this->memory_.x_0[level-1];
                coarse_reduction -= this->memory_.energy[level-1]; 

                this->transfer(level-1).interpolate(this->memory_.s[level-1], this->memory_.s[level]);

                if(!this->skip_BC_checks()){
                    this->zero_correction_related_to_equality_constrain(this->function(level), this->memory_.s[level]);
                }

                E_old = this->memory_.energy[level]; 
                this->memory_.x[level] += this->memory_.s[level];

                this->compute_s_global(level, this->memory_.s_working[level]);
                E_new = this->get_multilevel_energy(this->function(level), this->memory_.s_working[level], level);

                //----------------------------------------------------------------------------
                //                        trial point acceptance
                //----------------------------------------------------------------------------
                ared = E_old - E_new;
                rho = ared / coarse_reduction;
                if(coarse_reduction<=0){
                    rho = 0;
                }

                bool coarse_corr_taken = false;
                if(rho > this->rho_tol())
                {
                    coarse_corr_taken = true;
                    this->memory_.energy[level] = E_new; 
                }
                else
                {
                    this->memory_.x[level] -= this->memory_.s[level];
                    this->compute_s_global(level, this->memory_.s_working[level]);
                }

                //----------------------------------------------------------------------------
                //                                  trust region update
                //----------------------------------------------------------------------------
                converged = this->delta_update(rho, level, this->memory_.s_working[level]);

                // because, x + Is_{l-1} does not have to be inside of the feasible set....
                // mostly case for rmtr_inf with bounds...
                if(rho > this->rho_tol() && converged==false){
                    converged = this->check_feasibility(level);
                }


                if(this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
                {
                    // just to see what is being printed
                    std::string status = "RMTR_coarse_corr_stat, level: " + std::to_string(level);
                    this->print_init_message(status, {" it. ", "   E_old     ", "   E_new", "ared   ",  "  coarse_level_reduction  ", "  rho  ", "  delta ", "taken"});
                    PrintInfo::print_iter_status(this->_it_global, {E_old, E_new, ared, coarse_reduction, rho, this->memory_.delta[level], Scalar(coarse_corr_taken) });
                }

                // terminate, since TR rad. does not allow to take more corrections on given level
                if(converged==true){
                    return true;
                }

            }
            else if(mpi_world_rank() ==0 && this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE)
            {
                std::cout<< "--------- Recursion terminated due to non-smoothness of the gradient ----------------------- \n";
            }

            //----------------------------------------------------------------------------
            //                        postsmoothing
            //----------------------------------------------------------------------------

            if(this->post_smoothing_steps()!=0)
            {
                // auto post_smoothing_solve_type = (!smoothness_flg) ? COARSE_SOLVE : POST_SMOOTHING;
                auto post_smoothing_solve_type = POST_SMOOTHING;
                this->local_tr_solve(level, post_smoothing_solve_type);
            }

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
        virtual bool local_tr_solve(const SizeType & level, const LocalSolveType & solve_type)
        {
            SizeType it_success = 0, it = 0;
            Scalar ared = 0. , pred = 0., rho = 0., energy_new=9e9;
            bool make_grad_updates = true, make_hess_updates = true, converged = false, delta_converged = false;

            const bool exact_solve_flg = (solve_type == COARSE_SOLVE) ? true : false;
            // std::cout<<"solve_type: "<< solve_type << "    exact_solve_flg: "<< exact_solve_flg << "  \n"; 

            this->initialize_local_solve(level, solve_type);


            // should be neccessary just first time, we enter given level
            // if(empty(this->memory_.s[level]))
            // {
            //     this->memory_.s[level] = local_zeros(local_size(this->memory_.x[level]));
            // }
            // TODO:: check 
            this->memory_.s[level].set(0.0);


            // TODO:: do check if this is only post-smoothing (check if s_working is initiailized otherwise to zero)
            // important, as this can be postsmoothing
            // also check if fine level 
            this->compute_s_global(level, this->memory_.s_working[level]);


            if(!(solve_type==PRE_SMOOTHING && level==this->n_levels()-1))
            {
                if( (solve_type==PRE_SMOOTHING && level < this->n_levels()-1) || (solve_type == COARSE_SOLVE))
                {
                    if(CONSISTENCY_LEVEL == FIRST_ORDER)
                    {
                        this->ml_derivs_.g[level] += this->ml_derivs_.g_diff[level]; 
                    }
                    else if(CONSISTENCY_LEVEL == GALERKIN)
                    {
                        this->ml_derivs_.g[level] = this->ml_derivs_.g_diff[level]; 
                    }
                    else if(CONSISTENCY_LEVEL == SECOND_ORDER)
                    {
                        this->ml_derivs_.g[level] += this->ml_derivs_.g_diff[level]; 

                        // memory_.H[level] = memory_.H[level] + memory_.H_diff[level]; 
                        this->ml_derivs_.H[level] += this->ml_derivs_.H_diff[level]; 
                        make_hess_updates = false;                
                    }                    
                    else
                    {
                        utopia_error("RMTR:: Consistency order not existent .... \n"); 
                    }
                }
                else
                {
                    this->get_multilevel_gradient(this->function(level), this->memory_.s_working[level], level);
                }

                // energy computations ... 
                if(solve_type != POST_SMOOTHING)
                {
                    if(CONSISTENCY_LEVEL == SECOND_ORDER)
                    {
                        this->function(level).value(this->memory_.x[level], this->memory_.energy[level]); 
                        this->memory_.energy[level] = this->get_multilevel_energy(this->function(level), this->memory_.s_working[level], level);
                    }
                    else
                    {
                        this->memory_.energy[level] = this->get_multilevel_energy(this->function(level), this->memory_.s_working[level], level);
                    }
                }

                this->memory_.gnorm[level] = this->criticality_measure(level);
            }

            converged  = this->check_local_convergence(it, it_success,  this->memory_.gnorm[level], level, this->memory_.delta[level], solve_type);

            if(this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
            {
                this->print_level_info(level);
                PrintInfo::print_iter_status(0, {this->memory_.gnorm[level],  this->memory_.energy[level], ared, pred, rho, this->memory_.delta[level] });
            }

            it++;

            while(!converged)
            {
                if(make_hess_updates)
                {
                    this->get_multilevel_hessian(this->function(level), level);
                }

            //----------------------------------------------------------------------------
            //     solving constrained system to get correction and  building trial point
            //----------------------------------------------------------------------------
                // obtain correction
                this->solve_qp_subproblem(level, exact_solve_flg);

                // predicted reduction based on model
                pred = this->get_pred(level);

                // building trial point
                this->memory_.x[level] += this->memory_.s[level];

                this->compute_s_global(level, this->memory_.s_working[level]);
                energy_new = this->get_multilevel_energy(this->function(level), this->memory_.s_working[level], level);
                ared =  this->memory_.energy[level] - energy_new;

                rho = (ared < 0) ? 0.0 : ared/pred;
                rho = (rho != rho) ? 0.0 : rho;


                // update in hessian approx ...
                // TODO:: could be done in more elegant way....
                this->update_level(level);

            //----------------------------------------------------------------------------
            //     acceptance of trial point
            //----------------------------------------------------------------------------
                // good reduction, accept trial point
                if (rho >= this->rho_tol())
                {
                    it_success++;
                    this->memory_.energy[level] = energy_new;
                    make_grad_updates = true;
                }
                else
                {
                    this->memory_.x[level] -= this->memory_.s[level]; // return iterate into its initial state
                    this->compute_s_global(level, this->memory_.s_working[level]);
                    make_grad_updates = false;
                }

                //----------------------------------------------------------------------------
                //     updating level (deltas, hessian approx - new vectors, ...)
                //----------------------------------------------------------------------------
                delta_converged = this->delta_update(rho, level, this->memory_.s_working[level]);

                if(this->norm_schedule()==OUTER_CYCLE && this->verbosity_level() < VERBOSITY_LEVEL_VERY_VERBOSE && (solve_type==POST_SMOOTHING || solve_type == COARSE_SOLVE) && check_iter_convergence(it, it_success, level, solve_type))
                {
                    make_grad_updates = false; 
                }

                // can be more efficient, see commented line below 
                make_hess_updates =   make_grad_updates; 

                if(make_grad_updates)
                {
                    // std::cout<<"grad updated... \n"; 
                    // Vector g_old = memory_.g[level];
                    this->get_multilevel_gradient(this->function(level), this->memory_.s_working[level], level);
                    this->memory_.gnorm[level] = this->criticality_measure(level);

                    // make_hess_updates =   this->update_hessian(memory_.g[level], g_old, s, H, rho, g_norm);
                }
                // else
                // {
                //     make_hess_updates = false;  
                // }

                if(this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
                {
                    PrintInfo::print_iter_status(it, {this->memory_.gnorm[level], this->memory_.energy[level], ared, pred, rho, this->memory_.delta[level]});
                }

                converged  = (delta_converged  == true) ? true : this->check_local_convergence(it, it_success,  this->memory_.gnorm[level], level, this->memory_.delta[level], solve_type);

                if(level == this->n_levels()-1){
                    converged  = (converged  == true || this->memory_.gnorm[level] < this->atol()) ? true : false;
                }

                it++;

            }

            if(this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0){
                std::cout<< this->def_;
            }


            bool level_quit = ((this->criticality_measure_termination(this->memory_.gnorm[level]) == true) || delta_converged) ? true : false;
            return level_quit;
        }



    protected:
        virtual bool check_feasibility(const SizeType & /*level */ )
        {
            return false;
        }


        virtual void init_memory(const SizeType & fine_local_size) override
        {
            std::cout<<"------------------- init RMTR l2 .......... \n"; 
            RMTRBase::init_memory(fine_local_size); 

            const std::vector<SizeType> & dofs =  this->local_level_dofs(); 
            
            // init deltas to some default value...
            for(Scalar l = 0; l < this->n_levels(); l ++){
                this->memory_.delta[l] = this->delta0();
                this->_tr_subproblems[l]->init_memory(dofs[l]); 
            }
        }


        virtual Scalar get_pred(const SizeType & level)
        {
            this->memory_.help[level] = this->ml_derivs_.H[level] * this->memory_.s[level]; 
            return (-1.0 * dot(this->ml_derivs_.g[level], this->memory_.s[level]) -0.5 *dot(this->memory_.help[level], this->memory_.s[level]));
        }


        virtual void init_level(const SizeType & level)
        {
            this->memory_.delta[level]  = this->memory_.delta[level+1];
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
                 intermediate_delta = std::max(this->gamma1() * this->memory_.delta[level], 1e-15);
            else if (rho > this->eta2() )
                 intermediate_delta = std::min(this->gamma2() * this->memory_.delta[level], 1e15);
            else
                intermediate_delta = this->memory_.delta[level];


            // on the finest level we work just with one radius
            if(level==this->n_levels()-1)
            {
                this->memory_.delta[level] = intermediate_delta;
                return false;
            }
            else
            {
                Scalar corr_norm = this->level_dependent_norm(s_global, level);
                bool converged = this->delta_termination(corr_norm, level+1);

                if(converged && this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0){
                    std::cout<<"termination  due to small radius on level: "<< level << ". \n";
                }

                corr_norm = this->memory_.delta[level+1] - corr_norm;
                corr_norm = std::min(intermediate_delta, corr_norm);

                if(corr_norm <= 0.0)
                    corr_norm = 0.0;

                this->memory_.delta[level] = corr_norm;

                return converged;
            }
        }


        virtual bool update_level(const SizeType & /*level*/)
        {
            return false;
        }


        virtual void initialize_local_solve(const SizeType & /*level*/, const LocalSolveType & /*solve_type*/)
        {


        }


        /**
         * @brief      Computes norm of coarse level vector wrt to fine level
         *
         * @param[in]  u          The current iterate
         * @param[in]  current_l  The current level
         *
         */
        virtual Scalar level_dependent_norm(const Vector & u, const SizeType & level)
        {
            if(level == this->n_levels()-1){
                return norm2(u);
            }
            else
            {
                this->transfer(level).interpolate(u, this->memory_.help[level+1]);
                return norm2(this->memory_.help[level+1]);
            }
        }

        // ---------------------------------- convergence checks -------------------------------
        virtual bool check_global_convergence(const SizeType & it, const Scalar & r_norm, const Scalar & rel_norm, const Scalar & delta) override
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
        virtual bool check_local_convergence(const SizeType & it, const SizeType & it_success, const Scalar & g_norm, const SizeType & level, const Scalar & delta, const LocalSolveType & solve_type)
        {
            if(check_iter_convergence(it, it_success, level, solve_type))
                return true; 
            else if(delta < this->delta_min())
                return true;

            return this->criticality_measure_termination(g_norm);
        }


        /**
         * @brief      Checks for termination due to max it 
         *
         * @param[in]  it_success  Number of succesfull iterations on given level
         * @param[in]  g_norm      Norm of gradient
         * @param[in]  level       The level
         * @param[in]  delta       The delta
         *
         */
        virtual bool check_iter_convergence(const SizeType & it, const SizeType & it_success, const SizeType & level, const LocalSolveType & solve_type)
        {
            // coarse one
            if(level == 0 && (it_success >= this->max_sucessful_coarse_it() || it >= this->max_coarse_it()))
            {
                return true;
            }
            // every other level
            else if (level > 0 && solve_type == PRE_SMOOTHING)
            {
                if(it >= this->pre_smoothing_steps() || it_success >= this->max_sucessful_smoothing_it())
                {
                    return true;
                }
            }
            else if (level > 0 && solve_type == POST_SMOOTHING)
            {
                if(it >= this->post_smoothing_steps() || it_success >= this->max_sucessful_smoothing_it())
                {
                    return true;
                }
            }

            return false;
        }


        /**
         * @brief      THis check guarantees that iterates at a lower level remain in the TR radius defined at the finer level
         *
         * @param[in]  corr_norm  The norm of sum of all corrections on given level
         * @param[in]  level      The level
         *
         */
        virtual bool delta_termination(const Scalar & corr_norm, const SizeType & level)
        {
            return (corr_norm > (1.0 - this->eps_delta_termination()) * this->memory_.delta[level]) ? true : false;
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
            return (g_norm < this->atol()) ? true : false;
        }


        virtual Scalar criticality_measure(const SizeType & level) override
        {
            return norm2(this->ml_derivs_.g[level]);
        }



        /**
         * @brief      "Heuristics", which decides if it makes sense to go to the coarse level or no
         *
         * @param[in]  g_restricted  Restricted gradient
         * @param[in]  g_coarse      Coarse level gradient
         *
         */
        virtual bool grad_smoothess_termination(const Vector & g_restricted, const Vector & g_coarse, const SizeType & /*level*/)
        {
            Scalar Rg_norm, g_norm;
            norms2(g_restricted, g_coarse, Rg_norm, g_norm);
            return (Rg_norm >= this->get_grad_smoothess_termination() * g_norm) ? true : false;
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
            if(rho > 0 && rho < this->hessian_update_eta())
                return true;

            Vector help = g_new - g_old - H * s;

            // Hessian approx is relativelly poor
            return (norm2(help) > this->hessian_update_delta() * g_norm) ? true : false;
        }


//----------------------------- QP solve -----------------------------------------------------------------

        /**
         * @brief      Solves TR subroblem for given level
         *
         * @param[in]  level  The level
         * @param[in]  flg  The exact solve flag
         *
         */
        virtual bool solve_qp_subproblem(const SizeType & level, const bool & flg)
        {
            // this params should not be as hardcodded as they are...
            _tr_subproblems[level]->atol(1e-14);
            
            if(flg){
                _tr_subproblems[level]->max_it(this->max_QP_coarse_it());
            }
            else{
                _tr_subproblems[level]->max_it(this->max_QP_smoothing_it());
            }

            this->memory_.s[level].set(0.0); 
            _tr_subproblems[level]->current_radius(this->memory_.delta[level]);
            
            // maybe create help vector for this
            this->ml_derivs_.g[level] *= - 1.0; 
            _tr_subproblems[level]->solve(this->ml_derivs_.H[level], this->ml_derivs_.g[level], this->memory_.s[level]);
            this->ml_derivs_.g[level] *= - 1.0; 
            
            return true;
        }


        virtual void compute_s_global(const SizeType & level, Vector & s_global)
        {
            if(empty(this->memory_.x_0[level]))
            {
                s_global.set(0.0);
            }
            else if(level < this->n_levels()-1)
            {
                s_global = this->memory_.x[level] - this->memory_.x_0[level];
            }
        }


    private:
        std::vector<TRSubproblemPtr>        _tr_subproblems;

    };

}

#endif //UTOPIA_RMTR_HPP