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

    enum LocalSolveType{    PRE_SMOOTHING  = 1,
                            POST_SMOOTHING = 2,
                            COARSE_SOLVE   = 0};



    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
    class RMTR : public NonlinearMultiLevelBase<Matrix, Vector>,
                 public TrustRegionBase<Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        typedef utopia::TRSubproblem<Matrix, Vector>        TRSubproblem;
        typedef std::shared_ptr<TRSubproblem>               TRSubproblemPtr;

        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef utopia::Level<Matrix, Vector>               Level;

        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;

    public:

        using TrustRegionBase<Vector>::delta_update;
        using TrustRegionBase<Vector>::get_pred;

       /**
        * @brief      RMTR base-class
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level.
        */
        RMTR(   const SizeType & n_levels):
                NonlinearMultiLevelBase<Matrix,Vector>(n_levels),
                _it_global(0),
                _skip_BC_checks(false),
                _max_coarse_it(2),
                _max_smoothing_it(1),
                _max_QP_smoothing_it(5),
                _max_QP_coarse_it(50),
                _eps_delta_termination(0.001),
                _grad_smoothess_termination(0.5),
                _check_gradient_smoothness(true),
                _eps_grad_termination(1e-8),
                _hessian_update_delta(0.15),
                _hessian_update_eta(0.5),
                _verbosity_level(VERBOSITY_LEVEL_NORMAL),
                red_(FG_LIGHT_MAGENTA),
                def_(FG_DEFAULT),
                yellow_(FG_LIGHT_YELLOW),
                green_(FG_LIGHT_GREEN)
        {

        }


        virtual ~RMTR(){}

        virtual void read(Input &in) override
        {
            NonlinearMultiLevelBase<Matrix, Vector>::read(in);
            TrustRegionBase<Vector>::read(in);

            in.get("skip_BC_checks", _skip_BC_checks);
            in.get("max_coarse_it", _max_coarse_it);
            in.get("max_smoothing_it", _max_smoothing_it);
            in.get("max_QP_smoothing_it", _max_QP_smoothing_it);
            in.get("max_QP_coarse_it", _max_QP_coarse_it);

            in.get("eps_delta_termination", _eps_delta_termination);
            in.get("grad_smoothess_termination", _grad_smoothess_termination);
            in.get("check_gradient_smoothness", _check_gradient_smoothness);
            in.get("eps_grad_termination", _eps_grad_termination);

            in.get("hessian_update_delta", _hessian_update_delta);
            in.get("hessian_update_eta", _hessian_update_eta);
            // in.get("verbosity_level", _verbosity_level);

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
            NonlinearMultiLevelBase<Matrix, Vector>::print_usage(os);
            TrustRegionBase<Vector>::print_usage(os);

            this->print_param_usage(os, "skip_BC_checks", "bool", "Skip treatment of BC conditions.", "false");
            this->print_param_usage(os, "max_coarse_it", "int", "Maximum number of coarse iterations.", "2");
            this->print_param_usage(os, "max_smoothing_it", "int", "Maximum number of smoothing steps.", "1");
            this->print_param_usage(os, "max_QP_smoothing_it", "int", "Maximum number of iterations for fine level QP solver.", "5");
            this->print_param_usage(os, "max_QP_coarse_it", "int", "Maximum number of iterations for coarse level QP solver.", "50");

            this->print_param_usage(os, "eps_delta_termination", "real", "Constant used for quiting recursion based on size of tr. radius.", "0.001");
            this->print_param_usage(os, "grad_smoothess_termination", "real", "Constant used for quiting recursion based on smoothness of the gradient on given level.", "0.5");
            this->print_param_usage(os, "check_gradient_smoothness", "real", "Flag turning on/off check for gradient smoothiness.", "true");
            this->print_param_usage(os, "eps_grad_termination", "real", "Constant used for quiting recursion based on smoothness of the gradient on given level.", "1e-8");

            this->print_param_usage(os, "hessian_update_delta", "real", "Constant used for deciding whether to assemble fresh hessian or no.", "0.15");
            this->print_param_usage(os, "hessian_update_eta", "real", "Constant used for deciding whether to assemble fresh hessian or no.", "0.5");
            this->print_param_usage(os, "verbosity_level", "VerbosityLevel", "Specifies level of verbosity.", "VERBOSITY_LEVEL_NORMAL");

            this->print_param_usage(os, "coarse-QPSolver", "TRSubproblem", "Input parameters for fine level QP solvers.", "-");
            this->print_param_usage(os, "fine-QPSolver", "TRSubproblem", "Input parameters for coarse level QP solver.", "-");
        }


        virtual VerbosityLevel verbosity_level() const
        {
            return _verbosity_level;
        }

        virtual void verbosity_level(const VerbosityLevel & level )
        {

            _verbosity_level = this->verbose() ? level : VERBOSITY_LEVEL_QUIET;
        }


        virtual void set_grad_smoothess_termination(const Scalar & grad_smoothess_termination)
        {
            _grad_smoothess_termination = grad_smoothess_termination;
        }

        virtual Scalar  get_grad_smoothess_termination( ) const
        {
            return _grad_smoothess_termination;
        }


        virtual void skip_BC_checks(const bool skip_BC_checks)
        {
            _skip_BC_checks = skip_BC_checks;
        }

        virtual bool skip_BC_checks() const
        {
            return _skip_BC_checks;
        }

        using NonlinearMultiLevelBase<Matrix, Vector>::solve;

        virtual std::string name() override { return "RMTR";  }


        virtual void set_eps_grad_termination(const Scalar & eps_grad_termination)
        {
            _eps_grad_termination = eps_grad_termination;
        }


        virtual void max_coarse_it(const SizeType & max_coarse_it)
        {
            _max_coarse_it = max_coarse_it;
        }


        // this parameter define total number of smoothing its - sucessful and unsuccessful
        // while pre_smooting and post_smoothing count just sucessful iterations
        virtual void max_smoothing_it(const SizeType & max_smoothing_it)
        {
            _max_smoothing_it = max_smoothing_it;
        }

        virtual SizeType max_coarse_it() const
        {
            return _max_coarse_it;
        }


        virtual SizeType max_smoothing_it() const
        {
            return _max_smoothing_it;
        }


        virtual void max_QP_smoothing_it(const SizeType & num_it)
        {
            _max_QP_smoothing_it = num_it;
        }

        virtual void max_QP_coarse_it(const SizeType & num_it)
        {
            _max_QP_coarse_it = num_it;
        }

        virtual SizeType max_QP_coarse_it() const
        {
            return _max_QP_coarse_it;
        }

        virtual SizeType max_QP_smoothing_it() const
        {
            return _max_QP_smoothing_it;
        }


        virtual void check_grad_smoothness(const bool flg)
        {
            _check_gradient_smoothness = flg;
        }

        virtual bool check_grad_smoothness() const
        {
            return _check_gradient_smoothness;
        }


        virtual void handle_equality_constraints()
        {
            const SizeType L = this->n_levels();

            for(SizeType i = 0; i < L-1; ++i) {
                const auto &flags = this->function(i+1).get_eq_constrains_flg();
                assert(!empty(flags));

                this->transfer(i).handle_equality_constraints(flags);
            }
        }


        bool set_coarse_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy)
        {
            if(_tr_subproblems.size() != this->n_levels())
                _tr_subproblems.resize(this->n_levels());

            _tr_subproblems[0] = strategy;

            return true;
        }

        bool set_fine_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy)
        {
            if(_tr_subproblems.size() != this->n_levels())
                _tr_subproblems.resize(this->n_levels());

            // starting from level 1 ....
            for(std::size_t l = 1; l != this->n_levels(); ++l)
                _tr_subproblems[l] = std::shared_ptr<TRSubproblem>(strategy->clone());


            return true;
        }


        bool set_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy, const SizeType & level)
        {
            if(_tr_subproblems.size() != this->n_levels())
                _tr_subproblems.resize(this->n_levels());

            if(level <= this->n_levels())
            {
                _tr_subproblems[level] = strategy;
            }
            else{
                utopia_error("utopia::RMTR::set_tr_strategy:: Requested level exceeds number of levels in ML hierarchy. \n");
            }

            return true;
        }


        bool set_tr_strategies(const std::vector<TRSubproblemPtr> &strategies)
        {
            if(strategies.size() != this->n_levels()){
                utopia_error("utopia::RMTR::set_tr_strategies:: Number of tr strategies MUST be equal to number of levels in ML hierarchy. \n");
            }

            _tr_subproblems = strategies;

            return true;
        }


        /**
         * @brief      The solve function for multigrid method.
         *
         * @param[in]  rhs   The right hand side.
         * @param      x_0   The initial guess.
         *
         */
        virtual bool solve(Vector &x_h) override
        {
            if(!this->check_initialization())
                return false;


            bool converged = false;
            SizeType fine_level = this->n_levels()-1;
            Scalar r_norm, r0_norm, rel_norm;


            //-------------- INITIALIZATIONS ---------------
            if(!skip_BC_checks()){
                this->function(fine_level).get_eq_constrains_values(x_h);
            }

            SizeType fine_local_size = local_size(x_h).get(0);

            // this->status_.clear();
            this->init_memory(fine_local_size);

            // initialize();

            memory_.x[fine_level] = x_h;
            memory_.g[fine_level] = local_zeros(local_size(memory_.x[fine_level]));

            if(!skip_BC_checks()){
                this->make_iterate_feasible(this->function(fine_level), memory_.x[fine_level]);
            }

            this->function(fine_level).gradient(memory_.x[fine_level], memory_.g[fine_level]);
            this->function(fine_level).value(memory_.x[fine_level], memory_.energy[fine_level]);

            r0_norm = this->criticality_measure(fine_level);
            _it_global = 0;

            //----------------------------------------------

            if(verbosity_level() >= VERBOSITY_LEVEL_NORMAL && mpi_world_rank() == 0)
            {
                std::cout << red_;
                std::string name_id = this->name() + "     Number of levels: " + std::to_string(fine_level+1)  + "   \n Fine level local dofs: " + std::to_string(fine_local_size);
                this->init_solver(name_id, {" it. ", "|| g ||", "   E "});

                PrintInfo::print_iter_status(_it_global, {r0_norm, memory_.energy[fine_level]});
                std::cout << def_;
            }


            while(!converged)
            {
                if(this->cycle_type() == MULTIPLICATIVE_CYCLE)
                    this->multiplicative_cycle(fine_level);
                else{
                    std::cout<<"ERROR::UTOPIA_RMTR << unknown cycle type, solving in multiplicative manner ... \n";
                    this->multiplicative_cycle(fine_level);
                }

                #ifdef CHECK_NUM_PRECISION_mode
                    if(has_nan_or_inf(memory_.x[fine_level]) == 1)
                    {
                        memory_.x[fine_level] = local_zeros(local_size(memory_.x[fine_level]));
                        return true;
                    }
                #endif

                this->function(fine_level).gradient(memory_.x[fine_level], memory_.g[fine_level]);
                this->function(fine_level).value(memory_.x[fine_level], memory_.energy[fine_level]);

                r_norm = this->criticality_measure(fine_level);
                rel_norm = r_norm/r0_norm;

                _it_global++;

                if(this->verbose() && mpi_world_rank() == 0)
                {
                    std::cout << red_;
                    if(this->verbosity_level() > VERBOSITY_LEVEL_NORMAL){
                        this->print_init_message("RMTR OUTER SOLVE", {" it. ", "|| g ||", "|| g ||/|| g_0|" , "   E "});
                    }

                    PrintInfo::print_iter_status(_it_global, {r_norm, rel_norm,  memory_.energy[fine_level]});
                    std::cout << def_;
                }

                // check convergence
                converged = this->check_global_convergence(_it_global, r_norm, rel_norm, memory_.delta[fine_level]);
            }

            // benchmarking
            NonlinearMultiLevelBase<Matrix, Vector>::print_statistics(_it_global);
            x_h = memory_.x[fine_level];
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
        virtual bool multiplicative_cycle(const SizeType & level)
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
                converged = false;

            // making sure that correction does not exceed tr radius ...
            if(converged){
                return true;
            }

            this->compute_s_global(level, memory_.s_working[level]);
            this->get_multilevel_gradient(this->function(level), memory_.s_working[level], level);


            if(level == this->n_levels()-1)
            {
                converged =  this->criticality_measure_termination(this->criticality_measure(level));
                if(converged==true)
                    return true;
            }

            this->transfer(level-1).restrict(memory_.g[level], memory_.g_diff[level-1]);
            // std::cout<<"memory_.g_diff[level-1]: "<< norm2(memory_.g_diff[level-1]) << "  \n"; 



            this->transfer(level-1).project_down(memory_.x[level], memory_.x[level-1]);

            if(!skip_BC_checks()){
                this->make_iterate_feasible(this->function(level-1), memory_.x[level-1]);
            }

            //----------------------------------------------------------------------------
            //                  initializing coarse level (deltas, constraints, hessian approx, ...)
            //----------------------------------------------------------------------------
            this->init_level(level-1);

            //----------------------------------------------------------------------------
            //                   first order coarse level objective managment
            //----------------------------------------------------------------------------
            this->function(level-1).gradient(memory_.x[level-1], memory_.g[level-1]);
            // std::cout<<"memory_.g[level-1]: "<< norm2(memory_.g[level-1]) << "  \n"; 


            if(!skip_BC_checks())
            {
                if(CONSISTENCY_LEVEL != GALERKIN){
                    this->zero_correction_related_to_equality_constrain(this->function(level-1), memory_.g_diff[level-1]);
                }
            }

            if(this->check_grad_smoothness()){
                smoothness_flg = this->grad_smoothess_termination(memory_.g_diff[level-1], memory_.g[level-1], level-1);
            }
            else{
                smoothness_flg = true;
            }


            if(CONSISTENCY_LEVEL != GALERKIN)
            {
                memory_.g_diff[level-1] -= memory_.g[level-1];
                // std::cout<<"No-Galerkin:: memory_.g_diff[level-1]: "<< norm2(memory_.g_diff[level-1]) << "  \n"; 
            }

            //----------------------------------------------------------------------------
            //                   second order coarse level objective managment
            //----------------------------------------------------------------------------
            if(CONSISTENCY_LEVEL == SECOND_ORDER || CONSISTENCY_LEVEL == GALERKIN)
            {
                this->get_multilevel_hessian(this->function(level), level);
                this->transfer(level-1).restrict(memory_.H[level], memory_.H_diff[level-1]);

                if(CONSISTENCY_LEVEL == SECOND_ORDER)
                {
                    if(!skip_BC_checks()){
                        this->zero_correction_related_to_equality_constrain_mat(this->function(level-1), memory_.H_diff[level-1]);
                    }

                    this->function(level-1).hessian(memory_.x[level-1], memory_.H[level-1]);
                    memory_.H_diff[level-1] -=  memory_.H[level-1];
                }

            }

            //----------------------------------------------------------------------------
            //                   additional coarse level initialization...
            //----------------------------------------------------------------------------

            memory_.x_0[level-1]    = memory_.x[level-1];
            memory_.s[level-1]      = local_zeros(local_size(memory_.x[level-1]));

            // at this point s_global on coarse level is empty
            coarse_reduction = this->get_multilevel_energy(this->function(level-1), memory_.s[level-1], level-1);

            //----------------------------------------------------------------------------
            //               recursion  / Taylor correction
            //----------------------------------------------------------------------------
            if(level == 1 && smoothness_flg)
            {
                // std::cout<<"before calling recursion:: memory_.g_diff[level-1]: "<< norm2(memory_.g_diff[level-1]) << "  \n"; 
                // std::cout<<"before calling recursion:: memory_.g[level-1]: "<< norm2(memory_.g[level-1]) << "  \n"; 
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
                memory_.s[level-1] = memory_.x[level-1] - memory_.x_0[level-1];
                coarse_reduction -= this->get_multilevel_energy(this->function(level-1), memory_.s[level-1], level-1);

                this->transfer(level-1).interpolate(memory_.s[level-1], memory_.s[level]);

                if(!skip_BC_checks()){
                    this->zero_correction_related_to_equality_constrain(this->function(level), memory_.s[level]);
                }

                this->compute_s_global(level, memory_.s_working[level]);
                E_old = this->get_multilevel_energy(this->function(level), memory_.s_working[level], level);

                // new test for dbg mode
                memory_.x[level] += memory_.s[level];

                this->compute_s_global(level, memory_.s_working[level]);
                E_new = this->get_multilevel_energy(this->function(level), memory_.s_working[level], level);

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
                    memory_.x[level] -= memory_.s[level];
                    this->compute_s_global(level, memory_.s_working[level]);
                }

                //----------------------------------------------------------------------------
                //                                  trust region update
                //----------------------------------------------------------------------------
                converged = this->delta_update(rho, level, memory_.s_working[level]);

                // because, x + Is_{l-1} does not have to be inside of the feasible set....
                // mostly case for rmtr_inf with bounds...
                if(rho > this->rho_tol() && converged==false){
                    converged = this->check_feasibility(level);
                }


                // terminate, since TR rad. does not allow to take more corrections on given level
                if(converged==true){
                    return true;
                }


                if(this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
                {
                    // just to see what is being printed
                    std::string status = "RMTR_coarse_corr_stat, level: " + std::to_string(level);
                    this->print_init_message(status, {" it. ", "   E_old     ", "   E_new", "ared   ",  "  coarse_level_reduction  ", "  rho  ", "  delta ", "taken"});
                    PrintInfo::print_iter_status(_it_global, {E_old, E_new, ared, coarse_reduction, rho, memory_.delta[level], coarse_corr_taken });
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
            Scalar ared = 0. , pred = 0., rho = 0., energy_old=9e9, energy_new=9e9, g_norm=1.0;
            bool make_grad_updates = true, make_hess_updates = true, converged = false, delta_converged = false;

            const bool exact_solve_flg = (solve_type == COARSE_SOLVE) ? true : false;
            // std::cout<<"solve_type: "<< solve_type << "    exact_solve_flg: "<< exact_solve_flg << "  \n"; 


            this->initialize_local_solve(level, solve_type);


            // should be neccessary just first time, we enter given level
            if(empty(memory_.s[level]))
            {
                memory_.s[level] = local_zeros(local_size(memory_.x[level]));
            }

            // important, as this can be postsmoothing
            this->compute_s_global(level, memory_.s_working[level]);


            if(solve_type==PRE_SMOOTHING && level==this->n_levels()-1)
            {
                std::cout<<"skipping computation of gradient.... \n"; 
                energy_old = memory_.energy[level];
            }
            else
            {
                this->get_multilevel_gradient(this->function(level), memory_.s_working[level], level);
                energy_old = this->get_multilevel_energy(this->function(level), memory_.s_working[level], level);
            }

            
            g_norm = this->criticality_measure(level);

            converged  = this->check_local_convergence(it, it_success,  g_norm, level, memory_.delta[level], solve_type);

            if(this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
            {
                this->print_level_info(level);
                PrintInfo::print_iter_status(0, {g_norm, energy_old, ared, pred, rho, memory_.delta[level] });
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
                // correction needs to get prepared
                memory_.s[level] *= 0.0;
                this->solve_qp_subproblem(level, exact_solve_flg);

                // predicted reduction based on model
                pred = this->get_pred(level);

                // building trial point
                memory_.x[level] += memory_.s[level];

                this->compute_s_global(level, memory_.s_working[level]);
                energy_new = this->get_multilevel_energy(this->function(level), memory_.s_working[level], level);
                ared = energy_old - energy_new;

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
                    energy_old = energy_new;
                    make_grad_updates = true;
                }
                else
                {
                    memory_.x[level] -= memory_.s[level]; // return iterate into its initial state
                    this->compute_s_global(level, memory_.s_working[level]);
                    make_grad_updates = false;
                }

                //----------------------------------------------------------------------------
                //     updating level (deltas, hessian approx - new vectors, ...)
                //----------------------------------------------------------------------------
                delta_converged = this->delta_update(rho, level, memory_.s_working[level]);

                // can be more efficient, see below 
                make_hess_updates =   make_grad_updates; 

                if(make_grad_updates)
                {
                    Vector g_old = memory_.g[level];
                    this->get_multilevel_gradient(this->function(level), memory_.s_working[level], level);
                    g_norm = this->criticality_measure(level);

                    // make_hess_updates =   this->update_hessian(memory_.g[level], g_old, s, H, rho, g_norm);
                }

                if(this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
                    PrintInfo::print_iter_status(it, {g_norm, energy_new, ared, pred, rho, memory_.delta[level]});

                converged  = (delta_converged  == true) ? true : this->check_local_convergence(it, it_success,  g_norm, level, memory_.delta[level], solve_type);

                if(level == this->n_levels()-1){
                    converged  = (converged  == true || g_norm < this->atol()) ? true : false;
                }

                it++;

            }

            if(this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
                std::cout<< def_;


            bool level_quit = ((this->criticality_measure_termination(g_norm) == true) || delta_converged) ? true : false;
            return level_quit;
        }



    protected:

        virtual bool check_initialization()
        {
            if(this->level_functions_.size() != this->n_levels()){
                utopia_error("utopia::RMTR:: number of level Functions and levels not equal. \n");
                return false;
            }
            if(this->transfers_.size() + 1 != this->n_levels()){
                utopia_error("utopia::RMTR:: number of transfers and levels not equal. \n");
                return false;
            }

            // if(this->_tr_subproblems.size() != this->n_levels()){
            //     utopia_error("utopia::RMTR:: number of QP solvers and levels not equal. \n");
            //     return false;
            // }

            return true;
        }

        virtual bool check_feasibility(const SizeType & level )
        {
            return false;
        }

        virtual void init_memory(const SizeType & /*fine_local_size */) override
        {
            memory_.init(this->n_levels());

            // init deltas to some default value...
            for(Scalar l = 0; l < this->n_levels(); l ++)
                memory_.delta[l] = this->delta0();

        }


        virtual Scalar get_pred(const SizeType & level)
        {
            return (-1.0 * dot(memory_.g[level], memory_.s[level]) -0.5 *dot(memory_.H[level] * memory_.s[level], memory_.s[level]));
        }


        virtual void init_level(const SizeType & level)
        {
            memory_.delta[level]  = memory_.delta[level+1];
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
                 intermediate_delta = std::max(this->gamma1() * memory_.delta[level], 1e-15);
            else if (rho > this->eta2() )
                 intermediate_delta = std::min(this->gamma2() * memory_.delta[level], 1e15);
            else
                intermediate_delta = memory_.delta[level];


            // on the finest level we work just with one radius
            if(level==this->n_levels()-1)
            {
                memory_.delta[level] = intermediate_delta;
                return false;
            }
            else
            {
                Scalar corr_norm = this->level_dependent_norm(s_global, level);
                bool converged = this->delta_termination(corr_norm, level+1);

                if(converged && verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0){
                    std::cout<<"termination  due to small radius on level: "<< level+1 << ". \n";
                }

                corr_norm = memory_.delta[level+1] - corr_norm;
                corr_norm = std::min(intermediate_delta, corr_norm);

                if(corr_norm <= 0.0)
                    corr_norm = 0.0;

                memory_.delta[level] = corr_norm;

                return converged;
            }
        }


        virtual bool update_level(const SizeType & level)
        {
            return false;
        }


        virtual void initialize_local_solve(const SizeType & level, const LocalSolveType & solve_type)
        {


        }


        /**
         * @brief      Computes norm of coarse level vector wrt to fine level
         *
         * @param[in]  u          The current iterate
         * @param[in]  current_l  The current level
         *
         */
        virtual Scalar level_dependent_norm(const Vector & u, const SizeType & current_l)
        {
            if(current_l == this->n_levels()-1)
                return norm2(u);
            else
            {
                Vector s; // carries over prolongated correction
                this->transfer(current_l).interpolate(u, s);
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
        virtual bool check_local_convergence(const SizeType & it, const SizeType & it_success, const Scalar & g_norm, const SizeType & level, const Scalar & delta, const LocalSolveType & solve_type)
        {
            // coarse one
            if(level == 0 && (it_success >= _max_coarse_it))
            {
                return true;
            }
            // every other level
            else if (level > 0 && solve_type == PRE_SMOOTHING)
            {
                if(it_success >= this->pre_smoothing_steps() || it >= this->max_smoothing_it())
                    return true;
            }
            else if (level > 0 && solve_type == POST_SMOOTHING)
            {
                if(it_success >= this->post_smoothing_steps() || it >= this->max_smoothing_it())
                    return true;
            }
            else if(delta < this->delta_min())
                return true;

            return this->criticality_measure_termination(g_norm);
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
            return (corr_norm > (1.0 - _eps_delta_termination) * memory_.delta[level]) ? true : false;
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


        virtual Scalar criticality_measure(const SizeType & level)
        {
            return norm2(memory_.g[level]);
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
            _tr_subproblems[level]->atol(1e-16);
            
            if(flg){
                _tr_subproblems[level]->max_it(_max_QP_coarse_it);
            }
            else{
                _tr_subproblems[level]->max_it(_max_QP_smoothing_it);
            }
            // _tr_subproblems[level]->verbose("true"); 

            _tr_subproblems[level]->current_radius(memory_.delta[level]);
            _tr_subproblems[level]->solve(memory_.H[level], -1.0 * memory_.g[level], memory_.s[level]);

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
        virtual bool get_multilevel_hessian(const Fun & fun, const SizeType & level)
        {
            if(level < this->n_levels()-1){
                return MultilevelHessianEval<Matrix, Vector, CONSISTENCY_LEVEL>::compute_hessian(fun, memory_.x[level], memory_.H[level], memory_.H_diff[level]);
            }
            else{
                return fun.hessian(memory_.x[level], memory_.H[level]);
            }
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
        virtual bool get_multilevel_gradient(const Fun & fun, const Vector & s_global, const SizeType & level)
        {
            if(level < this->n_levels()-1)
            {
                return MultilevelGradientEval<Matrix, Vector, CONSISTENCY_LEVEL>::compute_gradient(fun, memory_.x[level], memory_.g[level], memory_.g_diff[level], memory_.H_diff[level], s_global);
            }
            else
            {
                // std::cout<<"--- fine level gradient eval.... \n";
                return fun.gradient(memory_.x[level], memory_.g[level]);
            }
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
        virtual Scalar get_multilevel_energy(const Fun & fun, const Vector & s_global, const SizeType & level)
        {
            if(level < this->n_levels()-1)
            {
                return MultilevelEnergyEval<Matrix, Vector, CONSISTENCY_LEVEL>::compute_energy(fun, memory_.x[level], memory_.g_diff[level], memory_.H_diff[level], s_global);
            }
            else
            {
                Scalar energy;
                fun.value(memory_.x[level], energy);
                return energy;
            }
        }


        virtual void compute_s_global(const SizeType & level, Vector & s_global)
        {
            if(empty( memory_.x_0[level]))
            {
                s_global = 0*memory_.x[level]; 
            }
            else if(level < this->n_levels()-1)
            {
                s_global = memory_.x[level] - memory_.x_0[level];
            }
        }


//----------------------------------------- printouts helpers -------------------------------------------

        /**
         * @brief      Prints some info related to level
         *
         * @param[in]  level  The level
         */
        virtual void print_level_info(const SizeType & level)
        {
            if(verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
            {
                if(level == 0)
                {
                    std::cout << yellow_;
                    std::string solver_type = "COARSE SOLVE:: " + std::to_string(level);
                    this->print_init_message(solver_type, {" it. ", "|| g ||", "   E + <g_diff, s>", "ared   ",  "  pred  ", "  rho  ", "  delta "});
                }
                else
                {
                    std::cout << green_;
                    std::string solver_type = "SMOOTHER:  " + std::to_string(level);
                    this->print_init_message(solver_type, {" it. ", "|| g ||", "   E + <g_diff, s>", "ared   ",  "  pred  ", "  rho  ", "  delta "});
                }
            }
        }

    private:
        std::vector<TRSubproblemPtr>        _tr_subproblems;

    protected:
        SizeType                            _it_global;                 /** * global iterate counter  */
        // ----------------------- PARAMETERS ----------------------
        bool                            _skip_BC_checks;

        SizeType                        _max_coarse_it;             /** * maximum iterations on coarse level   */
        SizeType                        _max_smoothing_it;          /** * max smoothing iterations  */

        SizeType                        _max_QP_smoothing_it;
        SizeType                        _max_QP_coarse_it;


        Scalar                         _eps_delta_termination;      /** * maximum delta allowed on coarse level - makes sure that coarse level corection stays inside fine level radius  */

        Scalar                         _grad_smoothess_termination; /** * determines when gradient is not smooth enough => does pay off to go to coarse level at all  */
        bool                           _check_gradient_smoothness;
        Scalar                         _eps_grad_termination;       /** * tolerance on grad  */

        Scalar                         _hessian_update_delta;       /** * tolerance used for updating hessians */
        Scalar                         _hessian_update_eta;         /** * tolerance used for updating hessians */

        VerbosityLevel                  _verbosity_level;

        ColorModifier red_;
        ColorModifier def_;
        ColorModifier yellow_;
        ColorModifier green_;


        RMTRLevelMemory <Matrix, Vector>         memory_;


    };

}

#endif //UTOPIA_RMTR_HPP