#ifndef UTOPIA_QUASI_RMTR_HPP
#define UTOPIA_QUASI_RMTR_HPP

#include "utopia_TRSubproblem.hpp"
#include "utopia_TRBase.hpp"

#include "utopia_MultiLevelEvaluations.hpp"
#include "utopia_RMTR.hpp"
#include "utopia_HessianApproximations.hpp"


namespace utopia 
{
    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
    class QuasiRMTR :   public RMTR<Matrix, Vector, CONSISTENCY_LEVEL>
    {
        typedef UTOPIA_SCALAR(Vector)      Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)   SizeType;

        typedef utopia::TRSubproblem<Matrix, Vector>                TRSubproblem; 
        typedef utopia::RMTR<Matrix, Vector, CONSISTENCY_LEVEL>     RMTR;

        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;

        typedef utopia::HessianApproximation<Vector>    HessianApproximation;
        typedef std::shared_ptr<HessianApproximation>   HessianApproxPtr; 


        static_assert(utopia::is_first_order<CONSISTENCY_LEVEL>::value, "utopia::QuasiRMTR does not support second order, nor galerkin consistency, nor Galerkin.");

    public:

        QuasiRMTR(  const SizeType & n_levels, 
                    const Parameters params = Parameters()): 
                    RMTR(n_levels, params)
        {
            set_parameters(params); 
            hessian_approxs_.resize(n_levels); 
        }

        virtual ~QuasiRMTR()
        {
            
        } 
        

        void set_parameters(const Parameters params) override
        {
            RMTR::set_parameters(params);    
        }


        virtual std::string name() override 
        { 
            return "Quasi_RMTR";  
        }
        

        virtual bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy)
        {
            if(hessian_approxs_.size() != this->n_levels())
                hessian_approxs_.resize(this->n_levels()); 

            for(auto l  = 0; l != hessian_approxs_.size(); ++l) 
                hessian_approxs_[l] = std::shared_ptr<HessianApproximation>(strategy->clone());

            return true;
        }
      
        virtual bool set_hessian_approximation_strategies(const std::vector<HessianApproxPtr> &strategies)
        {
            if(strategies.size() != this->n_levels()){
                utopia_error("utopia::QuasiRMTR::set_hessian_approximation_strategies:: Number of strategies does not equal with levels in ML hierarchy. \n"); 
            }
            
            hessian_approxs_ = strategies; 

            return true;
        }

        virtual bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy, const SizeType & level)
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


    protected:
        virtual void init_memory(const SizeType & fine_local_size) override 
        {   
            RMTR::init_memory(fine_local_size);
            hessian_approxs_[this->n_levels() - 1 ]->initialize();
        }

        virtual void init_level(const SizeType & level) override
        {
            RMTR::init_level(level);

            // maybe different strategies over here...
            // smoothed vectors and so on...
            hessian_approxs_[level]->initialize();
        }


        virtual bool get_multilevel_hessian(const Fun & fun, const SizeType & level) override
        {
            return false; 
        }


        virtual bool solve_qp_subproblem(const SizeType & level, const bool & flg) override
        {
            this->_tr_subproblems[level]->atol(1e-16);
            if(flg)
                this->_tr_subproblems[level]->max_it(this->_max_QP_coarse_it);
            else
                this->_tr_subproblems[level]->max_it(this->_max_QP_smoothing_it);


            auto multiplication_action = FunctionOperator<Vector>(hessian_approxs_[level]->get_apply_H()); 
            this->_tr_subproblems[level]->tr_constrained_solve(multiplication_action, this->memory_.g[level], this->memory_.s[level], this->memory_.delta[level]);            

            return true;
        }

        virtual Scalar get_pred(const SizeType & level) override
        {
            Scalar l_term = dot(this->memory_.g[level], this->memory_.s[level]);
            Scalar qp_term = hessian_approxs_[level]->compute_uHu_dot(this->memory_.s[level]); 
            return  (- l_term - 0.5 * qp_term); 
        }


        virtual bool update_level(const SizeType & level) override
        {
            Vector grad_old = this->memory_.g[level];
            this->get_multilevel_gradient(this->function(level), this->memory_.s_working[level], level);
            Vector y = this->memory_.g[level] - grad_old; 

            // swap back.... 
            this->memory_.g[level] = grad_old; 

            hessian_approxs_[level]->update(this->memory_.s[level], y);

            return true; 
        }


        virtual void reset_level(const SizeType & level) override
        {
            hessian_approxs_[level]->reset();
        }


        virtual bool check_initialization() override
        {
            RMTR::check_initialization();

            if(this->hessian_approxs_.size() != this->n_levels()){
                utopia_error("utopia::QuasiRMTR:: number of Hessian approximation streategies and levels not equal. \n"); 
                return false; 
            }

            return true; 
        }


    protected:   
        std::vector<HessianApproxPtr>  hessian_approxs_;


    };

}

#endif //UTOPIA_QUASI_RMTR_HPP