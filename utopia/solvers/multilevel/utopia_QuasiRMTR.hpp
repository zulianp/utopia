#ifndef UTOPIA_QUASI_RMTR_HPP
#define UTOPIA_QUASI_RMTR_HPP

#include "utopia_TRSubproblem.hpp"
#include "utopia_TRBase.hpp"

#include "utopia_MultiLevelEvaluations.hpp"
#include "utopia_RMTR.hpp"
#include "utopia_HessianApproximations.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class QuasiRMTR final:   public RMTRBase<Matrix, Vector, FIRST_ORDER_DF>
    {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        typedef utopia::RMTRBase<Matrix, Vector, FIRST_ORDER_DF>     RMTRBase;
        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun   Fun;

        using TRSubproblem = utopia::MatrixFreeTRSubproblem<Vector>;
        using TRSubproblemPtr = std::shared_ptr<TRSubproblem>;

        using HessianApproximation = utopia::HessianApproximation<Vector>;
        using HessianApproxPtr = std::shared_ptr<HessianApproximation>;

    public:
        QuasiRMTR(const SizeType & n_levels): RMTRBase(n_levels)
        {

        }

        ~QuasiRMTR() = default;

        void read(Input &in) override
        {
            RMTRBase::read(in);

            if(tr_subproblems_.size() > 0)
            {
                in.get("coarse-QPSolver", *tr_subproblems_[0]);

                for(auto i=1; i < tr_subproblems_.size(); i++)
                    in.get("fine-QPSolver", *tr_subproblems_[i]);
            }

            if(hessian_approxs_.size() > 0)
            {
                for(auto i=1; i < hessian_approxs_.size(); i++)
                    in.get("hessian-approx-strategy", *hessian_approxs_[i]);
            }
        }

        void print_usage(std::ostream &os) const override
        {
            RMTRBase::print_usage(os);

            this->print_param_usage(os, "coarse-QPSolver", "MatrixFreeTRSubproblem", "Input parameters for fine level QP solvers.", "-");
            this->print_param_usage(os, "fine-QPSolver", "MatrixFreeTRSubproblem", "Input parameters for coarse level QP solver.", "-");
            this->print_param_usage(os, "hessian-approx-strategy", "HessianApproximation", "Input parameters for hessian approximation strategies.", "-");
        }


        std::string name() override
        {
            return "Quasi_RMTR";
        }


        bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy)
        {
            if(hessian_approxs_.size() != this->n_levels())
                hessian_approxs_.resize(this->n_levels());

            for(auto l  = 0; l != hessian_approxs_.size(); ++l)
                hessian_approxs_[l] = std::shared_ptr<HessianApproximation>(strategy->clone());

            return true;
        }

        bool set_hessian_approximation_strategies(const std::vector<HessianApproxPtr> &strategies)
        {
            if(strategies.size() != this->n_levels()){
                utopia_error("utopia::QuasiRMTR::set_hessian_approximation_strategies:: Number of strategies does not equal with levels in ML hierarchy. \n");
            }

            hessian_approxs_ = strategies;

            return true;
        }

        bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy, const SizeType & level)
        {
            if(hessian_approxs_.size() != this->n_levels())
                hessian_approxs_.resize(this->n_levels());

            if(level <= this->n_levels())
            {
                hessian_approxs_[level] = strategy;
            }
            else{
                utopia_error("utopia::QuasiRMTR::set_tr_strategy:: Requested level exceeds number of levels in ML hierarchy. \n");
            }

            return true;
        }


        bool set_coarse_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy)
        {
            if(tr_subproblems_.size() != this->n_levels())
                tr_subproblems_.resize(this->n_levels());

            tr_subproblems_[0] = strategy;

            return true;
        }

        bool set_fine_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy)
        {
            if(tr_subproblems_.size() != this->n_levels())
                tr_subproblems_.resize(this->n_levels());

            // starting from level 1 ....
            for(std::size_t l = 1; l != tr_subproblems_.size(); ++l)
                tr_subproblems_[l] = std::shared_ptr<TRSubproblem>(strategy->clone());


            return true;
        }


        bool set_tr_strategies(const std::vector<TRSubproblemPtr> &strategies)
        {
            if(strategies.size() != this->n_levels()){
                utopia_error("utopia::RMTR::set_tr_strategies:: Number of tr strategies MUST be equal to number of levels in ML hierarchy. \n");
            }

            tr_subproblems_ = strategies;

            return true;
        }


    protected:
        void init_memory() override
        {
            const auto &layouts = this->local_level_layouts();
            bool same_fine_lo =  this->init_ ; 

            if(this->init_){
                same_fine_lo = layouts.back().same(layout(this->memory_.x.back())); 
            }

            if(!same_fine_lo){
                RMTRBase::init_memory();
                // const auto & layouts =  this->local_level_layouts();

                for(Scalar l = 0; l < this->n_levels(); l ++){
                    tr_subproblems_[l]->init_memory(layouts[l]);
                    hessian_approxs_[l]->initialize(this->memory_.x[l],this->ml_derivs_.g[l]);
                }
                this->init_ = true;   
            }
        }

        bool check_initialization() override
        {
            bool flg = RMTRBase::check_initialization();

            if(static_cast<SizeType>(tr_subproblems_.size()) != this->n_levels()){
                utopia_error("utopia::RMTR_l2_quasi:: number of level QP solvers and levels not equal. \n");
                flg = false;
            }

            if(static_cast<SizeType>(hessian_approxs_.size()) != this->n_levels()){
                utopia_error("utopia::RMTR_l2_quasi:: number of hessian approxiations and levels do not match. \n");
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
        bool delta_update(const Scalar & rho, const SizeType & level, const Vector & s_global) override
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



        bool solve_qp_subproblem(const SizeType & level, const bool & flg) override
        {
            Scalar atol_level = (level == this->n_levels()-1) ? this->atol() :  std::min(this->atol(), this->grad_smoothess_termination() * this->memory_.gnorm[level+1]);

            if (auto *tr_solver = dynamic_cast<IterativeSolver<Matrix, Vector> *>(tr_subproblems_[level].get())) {
                if(tr_solver->atol() > atol_level){
                    tr_solver->atol(atol_level);
                }

                if(flg){
                    tr_solver->max_it(this->max_QP_coarse_it());
                }
                else{
                    tr_solver->max_it(this->max_QP_smoothing_it());
                }
            }

            this->memory_.s[level].set(0.0);
            auto multiplication_action = hessian_approxs_[level]->build_apply_H();

            tr_subproblems_[level]->current_radius(this->memory_.delta[level]);
            this->ml_derivs_.g[level] *= - 1.0;
            tr_subproblems_[level]->solve(*multiplication_action, this->ml_derivs_.g[level], this->memory_.s[level]);
            this->ml_derivs_.g[level] *= - 1.0;

            if(has_nan_or_inf(this->memory_.s[level])){
                this->memory_.s[level].set(0.0);
            }

            return true;
        }


        // -------------------------- Matrix free stuff ----------------------------------------
        Scalar get_pred(const SizeType & level) override
        {
            UTOPIA_NO_ALLOC_BEGIN("RMTR::get_pred1");
            Scalar l_term = dot(this->ml_derivs_.g[level], this->memory_.s[level]);
            UTOPIA_NO_ALLOC_END();

            UTOPIA_NO_ALLOC_BEGIN("RMTR::get_pred2");
            Scalar qp_term = hessian_approxs_[level]->compute_uHu_dot(this->memory_.s[level]);
            UTOPIA_NO_ALLOC_END();

            return  (- l_term - 0.5 * qp_term);
        }


        bool update_level(const SizeType & level) override
        {
            this->memory_.help[level] = this->ml_derivs_.g[level];
            this->get_multilevel_gradient(this->function(level), level, this->memory_.s_working[level]);
            this->ml_derivs_.y[level] = this->ml_derivs_.g[level] - this->memory_.help[level];

            // swap back....
            this->ml_derivs_.g[level] = this->memory_.help[level];
            hessian_approxs_[level]->update(this->memory_.s[level], this->ml_derivs_.y[level], this->memory_.x[level], this->ml_derivs_.g[level]);

            return true;
        }

        void initialize_local_solve(const SizeType & level, const LocalSolveType & solve_type) override
        {
            // if(!(solve_type == PRE_SMOOTHING && level == this->n_levels()-1))
            // {
                if(solve_type == PRE_SMOOTHING || solve_type == COARSE_SOLVE)
                {
                    hessian_approxs_[level]->reset();
                    // hessian_approxs_[level]->initialize(this->memory_.x[level], this->ml_derivs_.g[level]);
                }
            // }
        }




    private:
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




    private:
        std::vector<HessianApproxPtr>       hessian_approxs_;
        std::vector<TRSubproblemPtr>        tr_subproblems_;


    };

}

#endif //UTOPIA_QUASI_RMTR_HPP