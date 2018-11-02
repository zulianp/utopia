#ifndef UTOPIA_QUASI_RMTR_HPP
#define UTOPIA_QUASI_RMTR_HPP

#include "utopia_TRSubproblem.hpp"
#include "utopia_TRBase.hpp"

#include "utopia_MultiLevelEvaluations.hpp"
#include "utopia_RMTR.hpp"
#include "utopia_HessianApprox.hpp"


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

        QuasiRMTR(  const std::shared_ptr<TRSubproblem> &tr_subproblem_coarse,  
                    const std::shared_ptr<TRSubproblem> &tr_subproblem_smoother,  
                    const Parameters params = Parameters()): 
                    RMTR(tr_subproblem_coarse, tr_subproblem_smoother, params)
        {
            set_parameters(params); 
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
            hessian_approx_clonable_ = strategy; 


            if(this->n_levels() > 0)
            {
                if(this->n_levels() != hessian_approxs_.size())
                    hessian_approxs_.resize(this->n_levels()); 

                for(auto l  = 0; l != hessian_approxs_.size(); ++l) 
                    hessian_approxs_[l] = std::shared_ptr<HessianApproximation>(strategy->clone());
            }
            else
            {
                utopia_error("utopia::QuasiRMTR::set_hessian_approximation_strategies:: Number of levels in ML hierarchy not set yet. \n Please set and call routine again .... \n"); 
            }

            return true;
        }
      
        virtual bool set_hessian_approximation_strategies(const std::vector<HessianApproxPtr> &strategies)
        {
            if(this->n_levels() > 0 && strategies.size() == this->n_levels())
            {
                if(this->n_levels() != hessian_approxs_.size())
                    hessian_approxs_.resize(this->n_levels()); 
            }
            else
            {
                utopia_error("utopia::QuasiRMTR::set_hessian_approximation_strategies:: Number of levels in ML hierarchy not set yet. \n Please set and call routine again .... \n"); 
            }
            
            hessian_approxs_ = strategies; 
            hessian_approx_clonable_ = hessian_approxs_[this->n_levels()-1]; 

            return true;
        }

        virtual bool set_hessian_approximation_strategy(const std::shared_ptr<HessianApproximation> &strategy, const SizeType & level)
        {
            if(this->n_levels() > 0 && hessian_approxs_.size() != this->n_levels())
            {
                hessian_approxs_.resize(this->n_levels()); 
            }
            else
            {
                utopia_error("utopia::RMTR::set_tr_strategy:: Number of levels in ML hierarchy not set yet. \n Please set and call routine again .... \n"); 
            }
            
            if(level <= this->n_levels())
            {
                hessian_approxs_[level] = strategy;
            }
            else
                utopia_error("utopia::RMTR::set_tr_strategy:: requested level exceeds number of levels in ML hierarchy.... \n"); 

            return true;
        }



    protected:
        virtual void init_memory(const SizeType & fine_local_size) override 
        {   
            RMTR::init_memory(fine_local_size);


            if(this->n_levels() != hessian_approxs_.size())
            {
                if(hessian_approx_clonable_ == nullptr)
                {
                    utopia_error("utopia::QuasiRMTR:: Hessian approximation strategy not provided..."); 
                }
                else
                {
                    hessian_approxs_.resize(this->n_levels());

                    for(std::size_t l = 0; l != hessian_approxs_.size(); ++l) 
                        hessian_approxs_[l] = std::shared_ptr<HessianApproximation>(hessian_approx_clonable_->clone());
                }
            }

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


    protected:   
        std::vector<HessianApproxPtr>  hessian_approxs_;
        std::shared_ptr<HessianApproximation> hessian_approx_clonable_;


    };

}

#endif //UTOPIA_QUASI_RMTR_HPP