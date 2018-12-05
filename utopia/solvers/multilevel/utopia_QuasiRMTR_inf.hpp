#ifndef UTOPIA_QUASI_RMTR_INF_HPP
#define UTOPIA_QUASI_RMTR_INF_HPP

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


namespace utopia 
{
    /**
     * @brief      The class for Quasi RMTR in infinity norm... 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
    class QuasiRMTR_inf :   public RMTR<Matrix, Vector, CONSISTENCY_LEVEL>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        typedef utopia::MatrixFreeQPSolver<Vector>  TRSubproblem;
        typedef std::shared_ptr<TRSubproblem>       TRSubproblemPtr; 

        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef utopia::Level<Matrix, Vector>               Level;


        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;


        typedef utopia::BoxConstraints<Vector>                          BoxConstraints;
        typedef utopia::RMTR<Matrix, Vector, CONSISTENCY_LEVEL>         RMTR;


        typedef utopia::HessianApproximation<Vector>    HessianApproximation;
        typedef std::shared_ptr<HessianApproximation>   HessianApproxPtr; 


        static_assert(utopia::is_first_order<CONSISTENCY_LEVEL>::value, "utopia::QuasiRMTR does not support second order, nor galerkin consistency, nor Galerkin.");        

    public:


        QuasiRMTR_inf(  const SizeType & n_levels,  
                        const Parameters params = Parameters()): 
                    RMTR(n_levels, params), 
                    has_box_constraints_(false)
        {
            set_parameters(params); 
            hessian_approxs_.resize(n_levels); 
        }

        virtual ~QuasiRMTR_inf()
        {
            // do we need to destroy some memory or no??? 
        } 
        

        void set_parameters(const Parameters params) override
        {
            RMTR::set_parameters(params);    
        }


        virtual std::string name() override 
        { 
            return "QuasiRMTR_inf";  
        }
        
        /**
         * @brief      Sets the box constraints.
         *
         * @param      box   The box
         *
         */
        virtual void set_box_constraints(BoxConstraints & box)
        {
          box_constraints_ = box; 
          has_box_constraints_ = true; 
        }

        /**
         * @brief      Gets the box constraints.
         *
         * @return     The box constraints.
         */
        virtual BoxConstraints & get_box_constraints()
        {
          return box_constraints_; 
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
        virtual void init_memory(const SizeType & fine_local_size) override 
        {   
            RMTR::init_memory(fine_local_size);
            
            constraints_memory_.init(this->n_levels()); 

            const SizeType fine_level = this->n_levels()-1; 
            const Scalar inf = std::numeric_limits<Scalar>::infinity();

            if(has_box_constraints_)
            {
                if(box_constraints_.has_upper_bound())
                    constraints_memory_.x_upper[fine_level] = *box_constraints_.upper_bound(); 
                else
                    constraints_memory_.x_upper[fine_level] = local_values(fine_local_size, inf); 

                if(box_constraints_.has_lower_bound())
                    constraints_memory_.x_lower[fine_level] = *box_constraints_.lower_bound(); 
                else
                    constraints_memory_.x_lower[fine_level] = local_values(fine_local_size, -1.0 * inf); 


                constraints_memory_.active_upper[fine_level] = constraints_memory_.x_upper[fine_level];
                constraints_memory_.active_lower[fine_level] = constraints_memory_.x_lower[fine_level];                
            }
            else
            {
                constraints_memory_.active_upper[fine_level] = local_values(fine_local_size, inf); 
                constraints_memory_.active_lower[fine_level] = local_values(fine_local_size, -1.0 * inf);
            }

            // inherited tr bound constraints... 
            constraints_memory_.tr_upper[fine_level] = local_values(fine_local_size, inf); 
            constraints_memory_.tr_lower[fine_level] = local_values(fine_local_size, -1.0 * inf);

            // precompute norms of prolongation operators needed for projections of constraints... 
            for(auto l = 0; l < this->n_levels() -1; l++)
                constraints_memory_.P_inf_norm[l] = this->transfer(l).interpolation_inf_norm(); 


            hessian_approxs_[this->n_levels() - 1 ]->initialize();
        }



        // since TR bounds are weak bounds...  
        virtual bool check_feasibility(const SizeType & level ) override
        {
            bool terminate = false; 

            {   
                Read<Vector> ru(constraints_memory_.tr_upper[level]); 
                Read<Vector> rl(constraints_memory_.tr_lower[level]); 
                Read<Vector> rx(this->memory_.x[level]); 

                Range r = range(constraints_memory_.tr_upper[level]);

                for(SizeType i = r.begin(); i != r.end(); ++i)
                {
                    Scalar xi = this->memory_.x[level].get(i); 
                    Scalar li = constraints_memory_.tr_lower[level].get(i); 
                    Scalar ui = constraints_memory_.tr_upper[level].get(i); 

                   if(xi < li || xi > ui)
                        terminate = true; 
                }
            }

            return terminate; 
        }



        // this routine is correct only under assumption, that P/R/I have only positive elements ... 
        virtual void init_level(const SizeType & level) override
        {
            RMTR::init_level(level);

            const SizeType finer_level = level+1; 

            //----------------------------------------------------------------------------
            //     soft projection of tr bounds
            //----------------------------------------------------------------------------
            Vector tr_fine_last_lower = this->memory_.x[finer_level] - local_values(local_size(this->memory_.x[finer_level]).get(0), this->memory_.delta[finer_level]);   
            {   
                ReadAndWrite<Vector> rv(tr_fine_last_lower); 
                Read<Vector> rl(constraints_memory_.tr_lower[finer_level]); 

                Range r = range(tr_fine_last_lower);

                for(SizeType i = r.begin(); i != r.end(); ++i) 
                    tr_fine_last_lower.set(i, std::max(constraints_memory_.tr_lower[finer_level].get(i), tr_fine_last_lower.get(i))); 
            }
            this->transfer(level).project_down(tr_fine_last_lower, constraints_memory_.tr_lower[level]);



            Vector tr_fine_last_upper = this->memory_.x[finer_level] + local_values(local_size(this->memory_.x[finer_level]).get(0), this->memory_.delta[finer_level]);   
            {   
                ReadAndWrite<Vector> rv(tr_fine_last_upper); 
                Read<Vector> rl(constraints_memory_.tr_upper[finer_level]); 

                Range r = range(tr_fine_last_upper);

                for(SizeType i = r.begin(); i != r.end(); ++i) 
                    tr_fine_last_upper.set(i, std::min(constraints_memory_.tr_upper[finer_level].get(i), tr_fine_last_upper.get(i))); 
            }
            this->transfer(level).project_down(tr_fine_last_upper, constraints_memory_.tr_upper[level]);


            if(has_box_constraints_)
            {
                //----------------------------------------------------------------------------
                //     projection of variable bounds to the coarse level
                //----------------------------------------------------------------------------
                Vector lx =  (constraints_memory_.x_lower[finer_level] - this->memory_.x[finer_level]); 
                Scalar lower_multiplier = 1.0/constraints_memory_.P_inf_norm[level] * max(lx); 
                constraints_memory_.x_lower[level] = this->memory_.x[level] + local_values(local_size(this->memory_.x[level]).get(0), lower_multiplier);

                Vector ux =  (constraints_memory_.x_upper[finer_level] - this->memory_.x[finer_level]); 
                Scalar upper_multiplier = 1.0/constraints_memory_.P_inf_norm[level] * min(ux); 
                constraints_memory_.x_upper[level] = this->memory_.x[level] + local_values(local_size(this->memory_.x[level]).get(0), upper_multiplier);

                //----------------------------------------------------------------------------
                //     intersect bounds on the coarse level
                //----------------------------------------------------------------------------
                constraints_memory_.active_upper[level] = local_zeros(local_size(this->memory_.x[level]).get(0)); 
                constraints_memory_.active_lower[level] = local_zeros(local_size(this->memory_.x[level]).get(0)); 
                {   
                    Write<Vector>   rv(constraints_memory_.active_upper[level]), rw(constraints_memory_.active_lower[level]); 
                    Read<Vector>    rl(this->memory_.x[level]), rq(constraints_memory_.x_lower[level]), re(constraints_memory_.x_upper[level]), rr(constraints_memory_.tr_lower[level]), rt(constraints_memory_.tr_upper[level]); 

                    Range r = range(this->memory_.x[level]);

                    for(SizeType i = r.begin(); i != r.end(); ++i)
                    {
                        constraints_memory_.active_upper[level].set(i, std::min(constraints_memory_.tr_upper[level].get(i), constraints_memory_.x_upper[level].get(i))); 
                        constraints_memory_.active_lower[level].set(i, std::max(constraints_memory_.tr_lower[level].get(i), constraints_memory_.x_lower[level].get(i))); 
                    }
                }
            }
            else
            {
                constraints_memory_.active_upper[level] = constraints_memory_.tr_upper[level]; 
                constraints_memory_.active_lower[level] = constraints_memory_.tr_lower[level]; 
            }


            // maybe different strategies over here...
            // smoothed vectors and so on...
            hessian_approxs_[level]->initialize();

        }

        virtual bool get_multilevel_hessian(const Fun & fun, const SizeType & level) override
        {
            return false; 
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
        virtual bool delta_update(const Scalar & rho, const SizeType & level, const Vector & s_global) override
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
            get_projection(x_g, constraints_memory_.active_lower[level], constraints_memory_.active_upper[level], Pc); 
            Pc -= this->memory_.x[level]; 
            Scalar Rg_norm =  norm2(Pc); 


            x_g = this->memory_.x[level] - g_coarse; 
            get_projection(x_g, constraints_memory_.active_lower[level], constraints_memory_.active_upper[level], Pc); 
            Pc -= this->memory_.x[level]; 
            Scalar  g_norm =  norm2(Pc); 

            return (Rg_norm >= this->get_grad_smoothess_termination() * g_norm) ? true : false;   
        }


        // measuring wrt to feasible set... 
        virtual Scalar criticality_measure(const SizeType & level) override
        {        
            Vector Pc; 
            Vector x_g = this->memory_.x[level] - this->memory_.g[level]; 

            get_projection(x_g, constraints_memory_.active_lower[level], constraints_memory_.active_upper[level], Pc); 
            
            Pc -= this->memory_.x[level]; 
            return norm2(Pc); 
        }


        /**
         * @brief      Projection onto feasible set
         *
         * @param[in]  x    Iterate
         * @param[in]  lb   lower bound
         * @param[in]  ub   upper bound
         * @param[in]  Pc   projection
         *
         */
        void get_projection(const Vector & x, const Vector &lb, const Vector &ub, Vector & Pc)
        {
            Pc = local_zeros(local_size(x).get(0)); 
            {
                Read<Vector> r_ub(ub), r_lb(lb), r_x(x);
                Write<Vector> wv(Pc); 

                each_write(Pc, [ub, lb, x](const SizeType i) -> double { 
                          Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);  
                          if(li >= xi)
                            return li; 
                          else
                            return (ui <= xi) ? ui : xi; }   );
            }
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
            Vector l = constraints_memory_.active_lower[level] - this->memory_.x[level]; 
            {   
                ReadAndWrite<Vector> rv(l); 
                each_write(l, [radius, l](const SizeType i) -> double { 
                    return  (l.get(i) >= -1*radius)  ? l.get(i) : -1 * radius;  }   );
            }

            Vector u =  constraints_memory_.active_upper[level] - this->memory_.x[level]; 
            {   
                ReadAndWrite<Vector> rv(u); 
                  each_write(u, [radius, u](const SizeType i) -> double { 
                      return  (u.get(i) <= radius)  ? u.get(i) : radius; }   );
            }

            // generating constraints to go for QP solve 
            auto box = make_box_constaints(std::make_shared<Vector>(l), std::make_shared<Vector>(u));


            // // setting should be really parameters from outside ... 
            // this->_tr_subproblems[level]->atol(1e-16); 

            // if(flg)
            //     this->_tr_subproblems[level]->max_it(this->max_QP_coarse_it()); 
            // else
            //     this->_tr_subproblems[level]->max_it(this->max_QP_smoothing_it());


            auto multiplication_action = hessian_approxs_[level]->build_apply_H(); 
            _tr_subproblems[level]->set_box_constraints(box); 

            this->_tr_subproblems[level]->solve(*multiplication_action, -1.0 * this->memory_.g[level], this->memory_.s[level]);


            // if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->_tr_subproblems[level].get()))
            //         tr_subproblem->tr_constrained_solve(*multiplication_action, this->memory_.g[level], this->memory_.s[level], box);            

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


    private:
        std::vector<TRSubproblemPtr>            _tr_subproblems; 

    protected:   
        ConstraintsLevelMemory <Vector>         constraints_memory_;

        BoxConstraints box_constraints_;        // constraints on the finest level.... 
        bool has_box_constraints_;               // as we can run rmtr with inf. norm also without constraints... 

        std::vector<HessianApproxPtr>  hessian_approxs_;

    };

}

#endif // UTOPIA_QUASI_RMTR_INF_HPP