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
    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
    class RMTR_l2 final: public RMTRBase<Matrix, Vector, CONSISTENCY_LEVEL>
    {
        typedef UTOPIA_SCALAR(Vector)           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)        SizeType;

        typedef utopia::TRSubproblem<Matrix, Vector>    TRSubproblem;
        typedef std::shared_ptr<TRSubproblem>           TRSubproblemPtr;

        typedef utopia::Transfer<Matrix, Vector>        Transfer;
        typedef utopia::Level<Matrix, Vector>           Level;

        typedef utopia::RMTRBase<Matrix, Vector, CONSISTENCY_LEVEL>     RMTRBase;
        typedef typename RMTRBase::Fun Fun;

        static_assert(!utopia::is_derivative_free<CONSISTENCY_LEVEL>::value, "utopia::RMTR_l2 does not support derivative free computations. Please use QuasiRMTR. ");

    public:

       /**
        * @brief      Implementation of RMTR with tr radius defined by l2 norm
        *
        * @param[in]  n_levels       Number of levels
        */
        RMTR_l2(const SizeType & n_levels): RMTRBase(n_levels)
        {

        }

         ~RMTR_l2(){}

        void read(Input &in) override
        {
            RMTRBase::read(in);

            const auto n_subproblems = _tr_subproblems.size();

            if(n_subproblems > 0)
            {
                in.get("coarse-QPSolver", *_tr_subproblems[0]);

                for(std::size_t i = 1; i < n_subproblems; i++)
                    in.get("fine-QPSolver", *_tr_subproblems[i]);
            }
        }

        void print_usage(std::ostream &os) const override
        {
            RMTRBase::print_usage(os);

            this->print_param_usage(os, "coarse-QPSolver", "TRSubproblem", "Input parameters for fine level QP solvers.", "-");
            this->print_param_usage(os, "fine-QPSolver", "TRSubproblem", "Input parameters for coarse level QP solver.", "-");
        }

        using NonlinearMultiLevelBase<Matrix, Vector>::solve;

        std::string name() override { return "RMTR"; }


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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    private:
        void init_memory() override
        {
            const auto &layouts = this->local_level_layouts();
            bool same_fine_lo =  this->init_ ; 

            if(this->init_){
                same_fine_lo = layouts.back().same(layout(this->memory_.x.back())); 
            }

            if(!same_fine_lo){
                RMTRBase::init_memory();
                // init deltas to some default value...
                for(Scalar l = 0; l < this->n_levels(); l ++){
                    this->_tr_subproblems[l]->init_memory(layouts[l]);
                }

                this->init_ = true; 
            }
        }

        bool check_initialization() override
        {
            bool flg = RMTRBase::check_initialization();

            if(static_cast<SizeType>(_tr_subproblems.size()) != this->n_levels()){
                utopia_error("utopia::RMTR_l2:: number of level QP solvers and levels not equal. \n");
                flg = false;
            }

            return flg;
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
        bool delta_update(const Scalar & rho, const SizeType & level, const Vector & s_global)
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

        /**
         * @brief      Computes norm of coarse level vector wrt to fine level
         *
         * @param[in]  u          The current iterate
         * @param[in]  current_l  The current level
         *
         */
        Scalar level_dependent_norm(const Vector & u, const SizeType & level)
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

        /**
         * @brief      THis check guarantees that iterates at a lower level remain in the TR radius defined at the finer level
         *
         * @param[in]  corr_norm  The norm of sum of all corrections on given level
         * @param[in]  level      The level
         *
         */
        bool delta_termination(const Scalar & corr_norm, const SizeType & level)
        {
            return (corr_norm > (1.0 - this->eps_delta_termination()) * this->memory_.delta[level]) ? true : false;
        }


        Scalar criticality_measure(const SizeType & level) override
        {
            return norm2(this->ml_derivs_.g[level]);
        }

        bool recursion_termination_smoothness(const Vector & g_restricted, const Vector & g_coarse, const SizeType & /*level*/) override
        {
            Scalar Rg_norm, g_norm;
            norms2(g_restricted, g_coarse, Rg_norm, g_norm);
            return (Rg_norm >= this->grad_smoothess_termination() * g_norm) ? true : false;
        }


//----------------------------- QP solve -----------------------------------------------------------------

        /**
         * @brief      Solves TR subroblem for given level
         *
         * @param[in]  level  The level
         * @param[in]  flg  The exact solve flag
         *
         */
        bool solve_qp_subproblem(const SizeType & level, const bool & flg) override
        {
            Scalar atol_level = (level == this->n_levels()-1) ? this->atol() :  std::min(this->atol(), this->grad_smoothess_termination() * this->memory_.gnorm[level+1]);
            if(_tr_subproblems[level]->atol() > atol_level){
                _tr_subproblems[level]->atol(atol_level);
            }

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

            if(has_nan_or_inf(this->memory_.s[level])){
                this->memory_.s[level].set(0.0);
            }

            return true;
        }

        Scalar get_pred(const SizeType & level) override
        {
            this->memory_.help[level] = this->ml_derivs_.H[level] * this->memory_.s[level];
            return (-1.0 * dot(this->ml_derivs_.g[level], this->memory_.s[level]) -0.5 *dot(this->memory_.help[level], this->memory_.s[level]));
        }

    private:
        std::vector<TRSubproblemPtr>        _tr_subproblems;

    };

}

#endif //UTOPIA_RMTR_HPP