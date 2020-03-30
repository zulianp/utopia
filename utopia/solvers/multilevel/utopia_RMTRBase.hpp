#ifndef UTOPIA_RMTR_BASE_HPP
#define UTOPIA_RMTR_BASE_HPP

#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"
#include "utopia_RMTRParams.hpp"

#include "utopia_TRSubproblem.hpp"
#include "utopia_Linear.hpp"
#include "utopia_Level.hpp"
#include "utopia_LS_Strategy.hpp"

#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_TRBase.hpp"

#include "utopia_FunEvalsIncludes.hpp"
#include "utopia_LevelMemory.hpp"

namespace utopia
{

    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
    class RMTRBase : public NonlinearMultiLevelBase<Matrix, Vector>, public RMTRParams<Vector>
    {
        public:
            typedef UTOPIA_SCALAR(Vector)                       Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;
            typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;

        RMTRBase(const SizeType & n_levels):    NonlinearMultiLevelBase<Matrix,Vector>(n_levels),
                                                RMTRParams<Vector>(),
                                                ml_derivs_(n_levels)
        {

        }

        virtual ~RMTRBase(){}

        virtual void read(Input &in) override
        {
            NonlinearMultiLevelBase<Matrix, Vector>::read(in);
            RMTRParams<Vector>::read(in);
        }

        virtual void print_usage(std::ostream &os) const override
        {
            NonlinearMultiLevelBase<Matrix, Vector>::print_usage(os);
            RMTRParams<Vector>::print_usage(os);
        }

    public:
        virtual bool solve(Vector &x_h) override;

    protected:
        virtual bool multiplicative_cycle(const SizeType & level);
        virtual bool local_tr_solve(const SizeType & level, const LocalSolveType & solve_type);


    protected:
        template<MultiLevelCoherence T = CONSISTENCY_LEVEL, enable_if_t<is_any<T, FIRST_ORDER, FIRST_ORDER_MGOPT, SECOND_ORDER, GALERKIN>::value, int> = 0 >
        bool get_multilevel_hessian(const Fun & fun, const SizeType & level)
        {
            return  ml_derivs_.compute_hessian(level, fun, memory_.x[level]);
        }

        template<MultiLevelCoherence T = CONSISTENCY_LEVEL, enable_if_t<is_any<T, FIRST_ORDER_DF>::value, int> = 0 >
        bool get_multilevel_hessian(const Fun & fun, const SizeType & level)
        {
            return  false;
        }

        virtual bool get_multilevel_gradient(const Fun & fun, const SizeType & level, const Vector & s_global) final{
            return ml_derivs_.compute_gradient(level, fun, memory_.x[level], s_global);
        }

        virtual bool get_multilevel_gradient(const Fun & fun, const SizeType & level) final{
            return ml_derivs_.compute_gradient(level, fun, memory_.x[level]);
        }


        virtual Scalar get_multilevel_energy(const Fun & fun, const SizeType & level, const Vector & s_global) final{
            return ml_derivs_.compute_energy(level, fun, memory_.x[level], s_global);
        }

        virtual Scalar get_multilevel_energy(const Fun & fun, const SizeType & level) final{
            return ml_derivs_.compute_energy(level, fun, memory_.x[level]);
        }


        virtual Scalar get_multilevel_gradient_energy(const Fun & fun, const SizeType & level, const Vector & s_global) final{
            return ml_derivs_.compute_gradient_energy(level, fun, memory_.x[level], s_global);
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        template<MultiLevelCoherence T = CONSISTENCY_LEVEL, enable_if_t<is_any<T, FIRST_ORDER, FIRST_ORDER_DF, FIRST_ORDER_MGOPT>::value, int> = 0 >
        bool init_consistency_terms(const SizeType & level)
        {
            //UTOPIA_NO_ALLOC_BEGIN("RMTR::region111");
            // Restricted fine level gradient
            this->transfer(level-1).restrict(this->ml_derivs_.g[level], this->ml_derivs_.g_diff[level-1]);

            // Projecting current iterate to obtain initial iterate on coarser grid
            this->transfer(level-1).project_down(this->memory_.x[level], this->memory_.x[level-1]);
            //UTOPIA_NO_ALLOC_END();

            //UTOPIA_NO_ALLOC_BEGIN("RMTR::region112");
            if(!this->skip_BC_checks()){
                this->make_iterate_feasible(this->function(level-1), this->memory_.x[level-1]);
            }
            //UTOPIA_NO_ALLOC_END();

            //----------------------------------------------------------------------------
            //    initializing coarse level (deltas, constraints, hessian approx, ...)
            //----------------------------------------------------------------------------
            UTOPIA_NO_ALLOC_BEGIN("RMTR::init_level");
            this->init_level(level-1);
            UTOPIA_NO_ALLOC_END();

            //----------------------------------------------------------------------------
            //                  first order coarse level objective managment
            //----------------------------------------------------------------------------
            //UTOPIA_NO_ALLOC_BEGIN("RMTR::region114");
            this->function(level-1).gradient(this->memory_.x[level-1], this->ml_derivs_.g[level-1]);
            //UTOPIA_NO_ALLOC_END();

            //UTOPIA_NO_ALLOC_BEGIN("RMTR::region1145");
            if(!this->skip_BC_checks()){
                this->zero_correction_related_to_equality_constrain(this->function(level-1), this->ml_derivs_.g_diff[level-1]);
            }
            //UTOPIA_NO_ALLOC_END();

            // UTOPIA_NO_ALLOC_BEGIN("RMTR::region115");
            bool smoothness_flg = this->check_grad_smoothness() ? this->recursion_termination_smoothness(this->ml_derivs_.g_diff[level-1], this->ml_derivs_.g[level-1], level-1) : true;
            this->ml_derivs_.g_diff[level-1] -= this->ml_derivs_.g[level-1];
            // UTOPIA_NO_ALLOC_END();

            return smoothness_flg;
        }

        template<MultiLevelCoherence T = CONSISTENCY_LEVEL, enable_if_t<is_same<T, SECOND_ORDER>::value, int> = 0 >
        bool init_consistency_terms(const SizeType & level)
        {
            UTOPIA_NO_ALLOC_BEGIN("RMTR::init_consistency_terms0");
            // Restricted fine level gradient
            this->transfer(level-1).restrict(this->ml_derivs_.g[level], this->ml_derivs_.g_diff[level-1]);

            // Projecting current iterate to obtain initial iterate on coarser grid
            this->transfer(level-1).project_down(this->memory_.x[level], this->memory_.x[level-1]);

            if(!this->skip_BC_checks()){
                this->make_iterate_feasible(this->function(level-1), this->memory_.x[level-1]);
            }

            //----------------------------------------------------------------------------
            //    initializing coarse level (deltas, constraints, hessian approx, ...)
            //----------------------------------------------------------------------------
            this->init_level(level-1);
            UTOPIA_NO_ALLOC_END();

            //----------------------------------------------------------------------------
            //                  first order coarse level objective managment
            //----------------------------------------------------------------------------
            UTOPIA_NO_ALLOC_BEGIN("RMTR::init_consistency_terms1");
            this->function(level-1).gradient(this->memory_.x[level-1], this->ml_derivs_.g[level-1]);

            if(!this->skip_BC_checks()){
                this->zero_correction_related_to_equality_constrain(this->function(level-1), this->ml_derivs_.g_diff[level-1]);
            }

            bool smoothness_flg = this->check_grad_smoothness() ? this->recursion_termination_smoothness(this->ml_derivs_.g_diff[level-1], this->ml_derivs_.g[level-1], level-1) : true;
            this->ml_derivs_.g_diff[level-1] -= this->ml_derivs_.g[level-1];
            UTOPIA_NO_ALLOC_END();

            //----------------------------------------------------------------------------
            //                   second order coarse level objective managment
            //----------------------------------------------------------------------------
            UTOPIA_NO_ALLOC_BEGIN("RMTR::hessian_comp2");
            this->get_multilevel_hessian(this->function(level), level);
            UTOPIA_NO_ALLOC_END();
            this->transfer(level-1).restrict(this->ml_derivs_.H[level], this->ml_derivs_.H_diff[level-1]);


            UTOPIA_NO_ALLOC_BEGIN("RMTR::init_consistency_terms22");
            if(!this->skip_BC_checks()){
                this->zero_correction_related_to_equality_constrain_mat(this->function(level-1), this->ml_derivs_.H_diff[level-1]);
            }
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("RMTR::hessian_comp3");
            this->function(level-1).hessian(this->memory_.x[level-1], this->ml_derivs_.H[level-1]);
            UTOPIA_NO_ALLOC_END();

            // memory_.H_diff[level-1] = memory_.H_diff[level-1] -  memory_.H[level-1];
            UTOPIA_NO_ALLOC_BEGIN("RMTR::init_consistency_terms24");
            this->ml_derivs_.H_diff[level-1] -= this->ml_derivs_.H[level-1];
            UTOPIA_NO_ALLOC_END();

            return smoothness_flg;

        }

        template<MultiLevelCoherence T = CONSISTENCY_LEVEL, enable_if_t<is_same<T, GALERKIN>::value, int> = 0 >
        bool init_consistency_terms(const SizeType & level)
        {
            UTOPIA_NO_ALLOC_BEGIN("RMTR::init_consistency_terms0");
            // Restricted fine level gradient
            this->transfer(level-1).restrict(this->ml_derivs_.g[level], this->ml_derivs_.g_diff[level-1]);

            if(!this->skip_BC_checks()){
                this->zero_correction_related_to_equality_constrain(this->function(level-1), this->ml_derivs_.g_diff[level-1]);
            }

            // Projecting current iterate to obtain initial iterate on coarser grid
            this->transfer(level-1).project_down(this->memory_.x[level], this->memory_.x[level-1]);

            //----------------------------------------------------------------------------
            //    initializing coarse level (deltas, constraints, hessian approx, ...)
            //----------------------------------------------------------------------------
            this->init_level(level-1);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("RMTR::hessian_comp2");
            this->get_multilevel_hessian(this->function(level), level);
            UTOPIA_NO_ALLOC_END();

            this->transfer(level-1).restrict(this->ml_derivs_.H[level], this->ml_derivs_.H_diff[level-1]);


            UTOPIA_NO_ALLOC_BEGIN("RMTR::init_consistency_terms1");
            if(!this->skip_BC_checks()){
                this->zero_correction_related_to_equality_constrain_mat(this->function(level-1), this->ml_derivs_.H_diff[level-1]);
            }
            UTOPIA_NO_ALLOC_END();

            return true;
        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        template<MultiLevelCoherence T = CONSISTENCY_LEVEL, enable_if_t<is_any<T, FIRST_ORDER, FIRST_ORDER_DF, FIRST_ORDER_MGOPT>::value, int> = 0 >
        bool  init_deriv_loc_solve(const Fun & fun, const SizeType & level, const LocalSolveType & solve_type)
        {
            if(!(solve_type==PRE_SMOOTHING && level==this->n_levels()-1)){

                if(solve_type==PRE_SMOOTHING || solve_type == COARSE_SOLVE){
                    this->ml_derivs_.g[level]  += this->ml_derivs_.g_diff[level];
                    this->memory_.gnorm[level] = this->criticality_measure(level);
                }
            }

            return true;
        }

        template<MultiLevelCoherence T = CONSISTENCY_LEVEL, enable_if_t<is_same<T, SECOND_ORDER>::value, int> = 0 >
        bool  init_deriv_loc_solve(const Fun & fun, const SizeType & level, const LocalSolveType & solve_type)
        {
            bool make_hess_updates = true;

            if(!(solve_type==PRE_SMOOTHING && level==this->n_levels()-1))
            {
                if( solve_type==PRE_SMOOTHING  || solve_type == COARSE_SOLVE)
                {
                    this->ml_derivs_.g[level] += this->ml_derivs_.g_diff[level];
                    this->memory_.gnorm[level] = this->criticality_measure(level);

                    this->ml_derivs_.H[level] += this->ml_derivs_.H_diff[level];
                    make_hess_updates = false;
                }
            }

            return make_hess_updates;
        }

        template<MultiLevelCoherence T = CONSISTENCY_LEVEL, enable_if_t<is_same<T, GALERKIN>::value, int> = 0 >
        bool  init_deriv_loc_solve(const Fun & fun, const SizeType & level, const LocalSolveType & solve_type)
        {
            if(!(solve_type==PRE_SMOOTHING && level==this->n_levels()-1)){

                if( (solve_type==PRE_SMOOTHING && level < this->n_levels()-1) || (solve_type == COARSE_SOLVE)){
                    this->ml_derivs_.g[level] = this->ml_derivs_.g_diff[level];
                    this->memory_.gnorm[level] = this->criticality_measure(level);
                }
            }

            return true;
        }

///////////////////////////////////////////////////////////// HELPERS ///////////////////////////////////////////////////////////////////////////////////////////////

        virtual void print_level_info(const SizeType & level) final
        {
            if(this->verbosity_level() >= VERBOSITY_LEVEL_VERY_VERBOSE && mpi_world_rank() == 0)
            {
                if(level == 0)
                {
                    std::cout << this->yellow_;
                    std::string solver_type = "COARSE SOLVE:: " + std::to_string(level);
                    this->print_init_message(solver_type, {" it. ", "|| g ||", "   E + <g_diff, x>", "ared   ",  "  pred  ", "  rho  ", "  delta "});
                }
                else
                {
                    std::cout << this->green_;
                    std::string solver_type = "SMOOTHER:  " + std::to_string(level);
                    this->print_init_message(solver_type, {" it. ", "|| g ||", "   E + <g_diff, x>", "ared   ",  "  pred  ", "  rho  ", "  delta "});
                }
            }
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
            if(rho > 0 && rho < this->hessian_update_eta())
                return true;

            // get rid of allocations
            Vector help = g_new - g_old - H * s;

            // Hessian approx is relativelly poor
            return (norm2(help) > this->hessian_update_delta() * g_norm) ? true : false;
        }


        virtual bool update_level(const SizeType & /*level*/)
        {
            return false;
        }


        virtual void compute_s_global(const SizeType & level, Vector & s_global) final
        {
            if(empty(this->memory_.x_0[level]))
            {
                utopia_error("this should not happen, remove when done testing... \n");
            }
            else if(level < this->n_levels()-1)
            {
                s_global = this->memory_.x[level] - this->memory_.x_0[level];
            }
            // we do not need to compute s_working on the fines level
        }

    ///////////////////////////////////////////// INITIALIZATIONS ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // to overwrite
        virtual void init_memory() override
        {
            const auto &layouts =  this->local_level_layouts();
            const std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector> > > & funs =  this->level_functions();

            ml_derivs_.init_memory(layouts, funs);
            memory_.init_memory(layouts);

            // init deltas to some default value...
            for(Scalar l = 0; l < this->n_levels(); l ++){
                this->memory_.delta[l] = this->delta0();
            }
        }

        void handle_equality_constraints()
        {
            const SizeType L = this->n_levels();

            for(SizeType i = 0; i < L-1; ++i) {
                const auto &flags = this->function(i+1).get_eq_constrains_flg();
                assert(!empty(flags));
                this->transfer(i).handle_equality_constraints(flags);
            }
        }

        virtual void init_level(const SizeType & level)
        {
            this->memory_.delta[level]  = this->memory_.delta[level+1];
        }

        virtual bool check_initialization()
        {
            if(static_cast<SizeType>(this->level_functions_.size()) != this->n_levels()){
                utopia_error("utopia::RMTR:: number of level Functions and levels not equal. \n");
                return false;
            }
            if(static_cast<SizeType>(this->transfers_.size()) + 1 != this->n_levels()){
                utopia_error("utopia::RMTR:: number of transfers and levels not equal. \n");
                return false;
            }

            return true;
        }

        virtual bool check_feasibility(const SizeType & /*level */ )
        {
            return false;
        }


    ///////////////////////////////////////////// CONVERGENCE CHECKS  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        bool check_global_convergence(const SizeType & it, const Scalar & r_norm, const Scalar & rel_norm, const Scalar & delta)
        {
            bool converged = NonlinearMultiLevelBase<Matrix, Vector>::check_convergence(it, r_norm, rel_norm, 1);

            if(delta < this->delta_min()){
                converged = true;
                this->exit_solver(it, ConvergenceReason::CONVERGED_TR_DELTA);
            }

            return converged;
        }


        virtual bool check_local_convergence(const SizeType & it, const SizeType & it_success, const SizeType & level, const Scalar & delta, const LocalSolveType & solve_type)
        {
            if(this->check_iter_convergence(it, it_success, level, solve_type)){
                return true;
            }
            else if(delta < this->delta_min()){
                return true;
            }

            return this->criticality_measure_termination(level);
        }


        virtual bool check_iter_convergence(const SizeType & it, const SizeType & it_success, const SizeType & level, const LocalSolveType & solve_type)
        {
            // coarse one
            if(level == 0 && (it_success >= this->max_sucessful_coarse_it() || it >= this->max_coarse_it())){
                return true;
            }
            // every other level
            else if (level > 0 && solve_type == PRE_SMOOTHING){
                if(it >= this->pre_smoothing_steps() || it_success >= this->max_sucessful_smoothing_it()){
                    return true;
                }
            }
            else if (level > 0 && solve_type == POST_SMOOTHING){
                if(it >= this->post_smoothing_steps() || it_success >= this->max_sucessful_smoothing_it()){
                    return true;
                }
            }

            return false;
        }


        virtual Scalar criticality_measure(const SizeType & level) = 0;
        virtual bool recursion_termination_smoothness(const Vector & g_restricted, const Vector & g_coarse, const SizeType & /*level*/) = 0;


        virtual bool criticality_measure_termination(const SizeType & level)
        {
                                                                                                                                       // this should be ideally different norm of grad on previous iteration
            Scalar atol_level = (level == this->n_levels()-1) ? this->atol() :  std::min(this->atol(), this->grad_smoothess_termination() * this->memory_.gnorm[level+1] );
            return (this->memory_.gnorm[level] < atol_level) ? true : false;
        }


///////////////////////////////////////////// TR-based stuff   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        virtual bool solve_qp_subproblem(const SizeType & level, const bool & flg) = 0;

        virtual void initialize_local_solve(const SizeType & /*level*/, const LocalSolveType & /*solve_type*/){ }

        virtual bool delta_update(const Scalar & rho, const SizeType & level, const Vector & s_global) = 0;

        virtual Scalar get_pred(const SizeType & level)  = 0;



    protected:
        RMTRLevelMemory<Matrix, Vector>         memory_;
        MultilevelDerivEval<Matrix, Vector, CONSISTENCY_LEVEL>  ml_derivs_;

    };

}

#endif //UTOPIA_RMTR_BASE_HPP