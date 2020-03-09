#ifndef UTOPIA_RMTR_INF_HPP
#define UTOPIA_RMTR_INF_HPP

#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"

#include "utopia_TRSubproblem.hpp"
#include "utopia_TrustRegionVariableBound.hpp"

#include "utopia_Linear.hpp"
#include "utopia_Level.hpp"

#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"
#include "utopia_TRBase.hpp"

#include "utopia_MultiLevelEvaluations.hpp"
#include "utopia_RMTR.hpp"

#include "utopia_IdentityTransfer.hpp"
#include "utopia_MultiLevelVariableBoundInterface.hpp"


namespace utopia
{
    
    /**
     * @brief      The class for RMTR in infinity norm...
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector, class MLConstraints,  MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
    class RMTR_inf final:   public RMTRBase<Matrix, Vector, CONSISTENCY_LEVEL>, public MLConstraints
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        typedef utopia::QPSolver<Matrix, Vector>            TRSubproblem;
        typedef std::shared_ptr<TRSubproblem>               TRSubproblemPtr;

        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef utopia::Level<Matrix, Vector>               Level;

        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;

        typedef utopia::BoxConstraints<Vector>                          BoxConstraints;
        typedef utopia::RMTRBase<Matrix, Vector, CONSISTENCY_LEVEL>         RMTR;


        // typedef MultilevelConstrInterface  MLConstraints;
        using MLConstraints::check_feasibility; 

    public:


       /**
        * @brief     RMTR with bound constraints ...
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level.
        */
        RMTR_inf(const SizeType & n_levels) : RMTRBase<Matrix, Vector, CONSISTENCY_LEVEL>(n_levels), MLConstraints(this->transfer())
                                                //has_box_constraints_(false) // not optional parameter
        {

        }

        void read(Input &in) override
        {
            RMTR::read(in);

            if(_tr_subproblems.size() > 0)
            {
                in.get("coarse-QPSolver", *_tr_subproblems[0]);

                for(auto i=1; i < static_cast<SizeType>(_tr_subproblems.size()); i++)
                    in.get("fine-QPSolver", *_tr_subproblems[i]);
            }
        }

        void print_usage(std::ostream &os) const override
        {
            RMTR::print_usage(os);

            this->print_param_usage(os, "coarse-QPSolver", "QPSolver", "Input parameters for fine level QP solvers.", "-");
            this->print_param_usage(os, "fine-QPSolver", "QPSolver", "Input parameters for coarse level QP solver.", "-");
        }

        ~RMTR_inf()
        {
            // do we need to destroy some memory or no???
        }


        std::string name() override
        {
            return "RMTR_inf";
        }

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
            for(auto l = 1; l != static_cast<SizeType>(_tr_subproblems.size()); ++l){
                _tr_subproblems[l] = std::shared_ptr<TRSubproblem>(strategy->clone());
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

    private:
        void init_memory() override
        {
            RMTR::init_memory();

            const std::vector<SizeType> & dofs =  this->local_level_dofs(); 
            MLConstraints::init_memory(dofs); 

            const SizeType fine_level = this->n_levels()-1;

            for(Scalar l = 0; l < this->n_levels(); l ++){
                _tr_subproblems[l]->init_memory(dofs[l]); 
            }               

            for(auto l = 0; l < fine_level; l++){
                this->transfer(l).init_memory(); 
            }
        }

        bool check_initialization() override
        {
            bool flg = RMTR::check_initialization(); 

            if(static_cast<SizeType>(_tr_subproblems.size()) != this->n_levels()){
                utopia_error("utopia::RMTR_inf:: number of level QP solvers and levels is not equal. \n");
                flg = false;
            }

            // if(static_cast<SizeType>(constraints_memory_.size()) != this->n_levels()){
            //     utopia_error("utopia::RMTR_l2_quasi:: number of hessian approxiations and levels do not match. \n");
            //     flg = false;
            // }            

            return flg; 
        }        


        Scalar get_pred(const SizeType & level) override
        {
            this->memory_.help[level] = this->ml_derivs_.H[level] * this->memory_.s[level]; 
            return (-1.0 * dot(this->ml_derivs_.g[level], this->memory_.s[level]) -0.5 *dot(this->memory_.help[level], this->memory_.s[level]));
        }


        bool check_feasibility(const SizeType & level ) override
        {
            return MLConstraints::check_feasibility(level, this->memory_.x[level]); 
        }


        // this routine is correct only under assumption, that P/R/I have only positive elements ...
        void init_level(const SizeType & level) override
        {
            RMTR::init_level(level);

            const SizeType finer_level = level+1;
            MLConstraints::init_level(level, this->memory_.x[finer_level], this->memory_.x[level], this->memory_.delta[finer_level]); 

            // let's see ... 
            // this->memory_.delta[level]  = this->delta0(); 

        }

        // -------------------------- tr radius managment ---------------------------------------------
        bool delta_update(const Scalar & rho, const SizeType & level, const Vector & /*s_global*/) override
        {
            Scalar intermediate_delta;

            // we could do also more sophisticated options, but lets not care for the moment ...
            if(rho < this->eta1()){
                intermediate_delta = std::max(this->gamma1() * this->memory_.delta[level], 1e-15);
            }
            else if (rho > this->eta2()){
                intermediate_delta = std::min(this->gamma2() * this->memory_.delta[level], 1e15);
            }
            else{
                intermediate_delta = this->memory_.delta[level];
            }

            this->memory_.delta[level] = intermediate_delta;

            return false;
        }


        bool recursion_termination_smoothness(const Vector & g_restricted, const Vector & g_coarse, const SizeType & level) override
        {
            // if we merge calls, reduction can be done together
            Scalar Rg_norm = MLConstraints::criticality_measure_inf(level, this->memory_.x[level], g_restricted); 
            Scalar g_norm = MLConstraints::criticality_measure_inf(level, this->memory_.x[level], g_coarse); 

            return (Rg_norm >= this->grad_smoothess_termination() * g_norm) ? true : false;
        }


        // measuring wrt to feasible set...
        Scalar criticality_measure(const SizeType & level) override
        {
            return MLConstraints::criticality_measure_inf(level, this->memory_.x[level], this->ml_derivs_.g[level]); 
        }

    public: // nvcc requires it to be public when using lambdas
        bool solve_qp_subproblem(const SizeType & level, const bool & flg) override
        {
            Scalar radius = this->memory_.delta[level];

            // first we need to prepare box of intersection of level constraints with tr. constraints
            std::shared_ptr<Vector> & lb = _tr_subproblems[level]->lower_bound(); 
            std::shared_ptr<Vector> & ub = _tr_subproblems[level]->upper_bound(); 

            const Vector & active_lower = this->active_lower(level); 
            const Vector & active_upper = this->active_upper(level); 

            // disp(active_lower, "active_lower"); 
            // disp(active_upper, "active_upper"); 

            *lb = active_lower - this->memory_.x[level];
            *ub = active_upper - this->memory_.x[level];
    
            {
                parallel_transform(*lb, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar {
                    return (xi >= -1.0*radius)  ? xi : -1.0*radius;
                });

                parallel_transform(*ub, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar{
                    return (xi <= radius)  ? xi : radius;
                });
            }            

            Scalar atol_level = (level == this->n_levels()-1) ? this->atol() :  std::min(this->atol(), this->grad_smoothess_termination() * this->memory_.gnorm[level+1]); 
            if(_tr_subproblems[level]->atol() > atol_level){
                _tr_subproblems[level]->atol(atol_level);  
            }

            if(flg){
                this->_tr_subproblems[level]->max_it(this->max_QP_coarse_it());
            }
            else{
                this->_tr_subproblems[level]->max_it(this->max_QP_smoothing_it());
            }

            this->ml_derivs_.g[level] *= - 1.0; 
            UTOPIA_NO_ALLOC_BEGIN("RMTR::qp_solve1");
            this->memory_.s[level].set(0.0); 
            this->_tr_subproblems[level]->solve(this->ml_derivs_.H[level], this->ml_derivs_.g[level], this->memory_.s[level]);
            UTOPIA_NO_ALLOC_END();
            this->ml_derivs_.g[level] *= - 1.0; 


            if(has_nan_or_inf(this->memory_.s[level])){
                this->memory_.s[level].set(0.0); 
            }
            else{
                // ----- just for debugging pourposes, to be commented out in the future... 
                MLConstraints::get_projection(*lb, *ub, this->memory_.s[level]); 
            }

            return true;
        }


    private:
        std::vector<TRSubproblemPtr>        _tr_subproblems;


    };

}

#endif //UTOPIA_RMTR_HPP
