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
    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
    class RMTR_inf :   public RMTR<Matrix, Vector, CONSISTENCY_LEVEL>, public MultilevelVariableBoundSolverInterface<Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        typedef utopia::QPSolver<Matrix, Vector>            TRSubproblem;
        typedef std::shared_ptr<TRSubproblem>               TRSubproblemPtr;

        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef utopia::Level<Matrix, Vector>               Level;

        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;

        typedef utopia::BoxConstraints<Vector>                          BoxConstraints;
        typedef utopia::RMTR<Matrix, Vector, CONSISTENCY_LEVEL>         RMTR;


        typedef utopia::MultilevelVariableBoundSolverInterface<Vector>  MLConstraints;
        using utopia::MultilevelVariableBoundSolverInterface<Vector>::check_feasibility; 

    public:


       /**
        * @brief     RMTR with bound constraints ...
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level.
        */
        RMTR_inf(   const SizeType & n_levels): RMTR(n_levels)//,
                                                //has_box_constraints_(false) // not optional parameter
        {

        }

        virtual void read(Input &in) override
        {
            RMTR::read(in);

            if(_tr_subproblems.size() > 0)
            {
                in.get("coarse-QPSolver", *_tr_subproblems[0]);

                for(auto i=1; i < static_cast<SizeType>(_tr_subproblems.size()); i++)
                    in.get("fine-QPSolver", *_tr_subproblems[i]);
            }
        }

        virtual void print_usage(std::ostream &os) const override
        {
            RMTR::print_usage(os);

            this->print_param_usage(os, "coarse-QPSolver", "QPSolver", "Input parameters for fine level QP solvers.", "-");
            this->print_param_usage(os, "fine-QPSolver", "QPSolver", "Input parameters for coarse level QP solver.", "-");
        }

        virtual ~RMTR_inf()
        {
            // do we need to destroy some memory or no???
        }


        virtual std::string name() override
        {
            return "RMTR_inf";
        }

        bool set_coarse_tr_strategy(const std::shared_ptr<TRSubproblem> &strategy)
        {
            if(static_cast<SizeType>(_tr_subproblems.size()) != this->n_levels())
                _tr_subproblems.resize(this->n_levels());

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



    protected:
        virtual void init_memory(const SizeType & fine_local_size) override
        {
            RMTR::init_memory(fine_local_size);
            MLConstraints::init_constr_memory(this->n_levels(), fine_local_size); 

            const SizeType fine_level = this->n_levels()-1;

            // precompute norms of prolongation operators needed for projections of constraints...
            for(auto l = 0; l < fine_level; l++){
                this->constraints_memory_.P_inf_norm[l] = this->transfer(l).interpolation_inf_norm();
            }

        }



        // since TR bounds are weak bounds...
        virtual bool check_feasibility(const SizeType & level ) override
        {
            return MLConstraints::check_feasibility(level, this->memory_.x[level]); 
        }



        // this routine is correct only under assumption, that P/R/I have only positive elements ...
        virtual void init_level(const SizeType & level) override
        {
            RMTR::init_level(level);

            const SizeType finer_level = level+1;
            MLConstraints::init_level(level, this->memory_.x[finer_level], this->memory_.x[level], this->memory_.delta[finer_level], this->transfer(level)); 

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
        virtual bool delta_update(const Scalar & rho, const SizeType & level, const Vector & /*s_global*/) override
        {
            Scalar intermediate_delta;

            // we could do also more sophisticated options, but lets not care for the moment ...
            if(rho < this->eta1())
                 intermediate_delta = std::max(this->gamma1() * this->memory_.delta[level], 1e-15);
            else if (rho > this->eta2() )
                 intermediate_delta = std::min(this->gamma2() * this->memory_.delta[level], 1e15);
            else
                intermediate_delta = this->memory_.delta[level];

            this->memory_.delta[level] = intermediate_delta;

            return false;
        }


        /**
         * @brief      "Heuristics", which decides if it makes sense to go to the coarse level or no
         *
         * @param[in]  g_restricted  Restricted gradient
         * @param[in]  g_coarse      Coarse level gradient
         *
         */
        virtual bool grad_smoothess_termination(const Vector & g_restricted, const Vector & g_coarse, const SizeType & level) override
        {
            Vector Pc;

            Vector x_g = this->memory_.x[level] - g_restricted;
            MLConstraints::get_projection(x_g, this->constraints_memory_.active_lower[level], this->constraints_memory_.active_upper[level], Pc);
            Pc -= this->memory_.x[level];
            Scalar Rg_norm =  norm2(Pc);


            x_g = this->memory_.x[level] - g_coarse;
            MLConstraints::get_projection(x_g, this->constraints_memory_.active_lower[level], this->constraints_memory_.active_upper[level], Pc);
            Pc -= this->memory_.x[level];
            Scalar  g_norm =  norm2(Pc);

            return (Rg_norm >= this->get_grad_smoothess_termination() * g_norm) ? true : false;
        }


        // measuring wrt to feasible set...
        virtual Scalar criticality_measure(const SizeType & level) override
        {
            return MLConstraints::criticality_measure_inf(level, this->memory_.x[level], this->ml_derivs_.g[level]); 
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
        virtual bool solve_qp_subproblem(const SizeType & level, const bool & flg) override
        {
            Scalar radius = this->memory_.delta[level];

            // first we need to prepare box of intersection of level constraints with tr. constraints
            Vector l = this->constraints_memory_.active_lower[level] - this->memory_.x[level];
            each_transform(l, l, [radius](const SizeType /*i*/, const Scalar val) -> Scalar {
                return (val >= -1*radius)  ? val : -1 * radius;  }
            );


            Vector u =  this->constraints_memory_.active_upper[level] - this->memory_.x[level];
            each_transform(u, u, [radius](const SizeType /*i*/, const Scalar val) -> Scalar {
              return (val <= radius)  ? val : radius; }
            );


            // generating constraints to go for QP solve
            auto box = make_box_constaints(std::make_shared<Vector>(l), std::make_shared<Vector>(u));


            // setting should be really parameters from outside ...
            this->_tr_subproblems[level]->atol(1e-14);

            if(flg){
                this->_tr_subproblems[level]->max_it(this->max_QP_coarse_it());
            }
            else{
                this->_tr_subproblems[level]->max_it(this->max_QP_smoothing_it());
            }


            _tr_subproblems[level]->set_box_constraints(box);
            this->_tr_subproblems[level]->solve(this->ml_derivs_.H[level], -1.0 * this->ml_derivs_.g[level], this->memory_.s[level]);



            // ----- just for debugging pourposes ----------------
            Vector s_old = this->memory_.s[level]; 
            MLConstraints::get_projection(s_old, l, u, this->memory_.s[level]); 


            return true;
        }


    private:
        std::vector<TRSubproblemPtr>        _tr_subproblems;


    };

}

#endif //UTOPIA_RMTR_HPP