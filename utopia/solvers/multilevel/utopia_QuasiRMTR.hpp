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

        typedef utopia::RMTR<Matrix, Vector, CONSISTENCY_LEVEL>         RMTR;
        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun   Fun;


        typedef utopia::MatrixFreeTRSubproblem<Vector>  TRSubproblem;
        typedef std::shared_ptr<TRSubproblem>           TRSubproblemPtr;


        typedef utopia::HessianApproximation<Vector>    HessianApproximation;
        typedef std::shared_ptr<HessianApproximation>   HessianApproxPtr;


        static_assert(utopia::is_first_order<CONSISTENCY_LEVEL>::value, "utopia::QuasiRMTR does not support second order, nor galerkin consistency, nor Galerkin.");

    public:

        QuasiRMTR(  const SizeType & n_levels):
                    RMTR(n_levels)
        {
         //   hessian_approxs_.resize(n_levels);
        }

        virtual ~QuasiRMTR()
        {

        }

        virtual void read(Input &in) override
        {
            RMTR::read(in);

            if(_tr_subproblems.size() > 0)
            {
                in.get("coarse-QPSolver", *_tr_subproblems[0]);

                for(auto i=1; i < _tr_subproblems.size(); i++)
                    in.get("fine-QPSolver", *_tr_subproblems[i]);
            }

            if(hessian_approxs_.size() > 0)
            {
                for(auto i=1; i < hessian_approxs_.size(); i++)
                    in.get("hessian-approx-strategy", *hessian_approxs_[i]);
            }
        }

        virtual void print_usage(std::ostream &os) const override
        {
            RMTR::print_usage(os);

            this->print_param_usage(os, "coarse-QPSolver", "MatrixFreeTRSubproblem", "Input parameters for fine level QP solvers.", "-");
            this->print_param_usage(os, "fine-QPSolver", "MatrixFreeTRSubproblem", "Input parameters for coarse level QP solver.", "-");
            this->print_param_usage(os, "hessian-approx-strategy", "HessianApproximation", "Input parameters for hessian approximation strategies.", "-");
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
            for(std::size_t l = 1; l != _tr_subproblems.size(); ++l)
                _tr_subproblems[l] = std::shared_ptr<TRSubproblem>(strategy->clone());


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


    protected:
        virtual void init_memory() override
        {
            RMTR::init_memory();

            SizeType fine_lev=this->n_levels() - 1; 
            hessian_approxs_[fine_lev]->initialize(this->memory_.x[fine_lev], this->memory_.g[fine_lev]);
        }


        virtual bool get_multilevel_hessian(const Fun & /*fun*/, const SizeType & /*level*/) override
        {
            return false;
        }


        virtual bool solve_qp_subproblem(const SizeType & level, const bool & flg) override
        {
            // TODO:: cast to iterative solvers
            // this->_tr_subproblems[level]->atol(1e-14);
            // if(flg){
            //     this->_tr_subproblems[level]->max_it(this->_max_QP_coarse_it);
            // }
            // else{
            //     this->_tr_subproblems[level]->max_it(this->_max_QP_smoothing_it);
            // }

            auto multiplication_action = hessian_approxs_[level]->build_apply_H();

            _tr_subproblems[level]->current_radius(this->memory_.delta[level]);
            _tr_subproblems[level]->solve(*multiplication_action, -1.0 * this->memory_.g[level], this->memory_.s[level]);

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

            hessian_approxs_[level]->update(this->memory_.s[level], y, this->memory_.x[level], this->memory_.g[level]);

            return true;
        }


        virtual void initialize_local_solve(const SizeType & level, const LocalSolveType & solve_type) override
        {
            if(!(solve_type == PRE_SMOOTHING && level == this->n_levels()-1))
            {
                // this is interesting heuristic
                if(solve_type == PRE_SMOOTHING || solve_type == COARSE_SOLVE)
                {
                    hessian_approxs_[level]->reset();
                    hessian_approxs_[level]->initialize(this->memory_.x[level], this->memory_.g[level]);                    
                }
            }
        }

        virtual Scalar criticality_measure(const SizeType & level) override
        {
            return norm2(this->ml_derivs_.g[level]);
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


    private:
        std::vector<TRSubproblemPtr>        _tr_subproblems;


    };

}

#endif //UTOPIA_QUASI_RMTR_HPP