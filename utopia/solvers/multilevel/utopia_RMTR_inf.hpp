#ifndef UTOPIA_RMTR_INF_HPP
#define UTOPIA_RMTR_INF_HPP

#include "utopia_NonLinearSmoother.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonlinearMultiLevelBase.hpp"

#include "utopia_TRSubproblem.hpp"
#include "utopia_TRBoxSubproblem.hpp"
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
     * @brief      The class for RMTR in infinity norm... 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector, MultiLevelCoherence CONSISTENCY_LEVEL = FIRST_ORDER>
    class RMTR_inf :   public RMTR<Matrix, Vector, CONSISTENCY_LEVEL>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        // pay attention that this one in inf norm... 
        typedef utopia::TRBoxSubproblem<Matrix, Vector>     TRSubproblem; 

        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef utopia::Level<Matrix, Vector>               Level;


        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;


        typedef utopia::BoxConstraints<Vector>                          BoxConstraints;
        typedef utopia::RMTR<Matrix, Vector, CONSISTENCY_LEVEL>         RMTR;

    public:


       /**
        * @brief     RMTR with bound constraints ... 
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level. 
        */
        RMTR_inf(   const std::shared_ptr<TRSubproblem> &tr_subproblem_coarse,  
                    const std::shared_ptr<TRSubproblem> &tr_subproblem_smoother,  
                    const Parameters params = Parameters()): 
                    
                    RMTR(tr_subproblem_coarse, tr_subproblem_smoother, params), 
                    has_box_constraints_(false)
        {
            set_parameters(params); 
        }

        virtual ~RMTR_inf()
        {
            // do we need to destroy some memory or no??? 
        } 
        

        void set_parameters(const Parameters params) override
        {
            RMTR::set_parameters(params);    
        }


        virtual std::string name() override 
        { 
            return "RMTR_inf";  
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


    protected:
        virtual void init_memory(const SizeType & fine_local_size) override 
        {   
            RMTR::init_memory(fine_local_size);
            constraints_memory_.init(this->n_levels()); 

            const SizeType fine_level = this->n_levels()-1; 
            const Scalar inf = std::numeric_limits<Scalar>::infinity();

            // bound constraints ....
            if(has_box_constraints_)
            {
                if(box_constraints_.has_upper_bound())
                    constraints_memory_.x_upper[fine_level] = *box_constraints_.upper_bound(); 
                else
                    constraints_memory_.x_upper[fine_level] = inf * local_values(fine_local_size, 1.0); 

                if(box_constraints_.has_lower_bound())
                    constraints_memory_.x_lower[fine_level] = *box_constraints_.lower_bound(); 
                else
                    constraints_memory_.x_lower[fine_level] = -1.0 * inf * local_values(fine_local_size, 1.0); 
            }
            else
            {
                constraints_memory_.x_upper[fine_level] = inf * local_values(fine_local_size, 1.0); 
                constraints_memory_.x_lower[fine_level] = -1.0 * inf * local_values(fine_local_size, 1.0);
            }   

            // inherited tr bound constraints... 
            constraints_memory_.tr_upper[fine_level] = inf * local_values(fine_local_size, 1.0); 
            constraints_memory_.tr_lower[fine_level] = -1.0 * inf * local_values(fine_local_size, 1.0);

            // active, for the finest level simply bound, since tr ones are artificial.... 
            constraints_memory_.active_upper[fine_level] = constraints_memory_.x_upper[fine_level];
            constraints_memory_.active_lower[fine_level] = constraints_memory_.x_lower[fine_level];


            // precompute norms of prolongation operators needed for projections of constraints... 
            for(auto l = 0; l < this->n_levels() -1; l++)
                constraints_memory_.P_inf_norm[l] = this->transfer(l).interpolation_inf_norm(); 

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
        virtual void init_coarse_level_constrains(const SizeType & level) override
        {
            // init delta on coarser level...
            this->memory_.delta[level-1]  = this->memory_.delta[level]; 

            //----------------------------------------------------------------------------
            //     soft projection of tr bounds
            //----------------------------------------------------------------------------
            Vector tr_fine_last_lower = this->memory_.x[level] - local_values(local_size(this->memory_.x[level]).get(0), this->memory_.delta[level]);   
            {   
                ReadAndWrite<Vector> rv(tr_fine_last_lower); 
                Read<Vector> rl(constraints_memory_.tr_lower[level]); 

                Range r = range(tr_fine_last_lower);

                for(SizeType i = r.begin(); i != r.end(); ++i) 
                    tr_fine_last_lower.set(i, std::max(constraints_memory_.tr_lower[level].get(i), tr_fine_last_lower.get(i))); 
            }
            this->transfer(level-1).project_down(tr_fine_last_lower, constraints_memory_.tr_lower[level-1]);



            Vector tr_fine_last_upper = this->memory_.x[level] + local_values(local_size(this->memory_.x[level]).get(0), this->memory_.delta[level]);   
            {   
                ReadAndWrite<Vector> rv(tr_fine_last_upper); 
                Read<Vector> rl(constraints_memory_.tr_upper[level]); 

                Range r = range(tr_fine_last_upper);

                for(SizeType i = r.begin(); i != r.end(); ++i) 
                    tr_fine_last_upper.set(i, std::min(constraints_memory_.tr_upper[level].get(i), tr_fine_last_upper.get(i))); 
            }
            this->transfer(level-1).project_down(tr_fine_last_upper, constraints_memory_.tr_upper[level-1]);


            //----------------------------------------------------------------------------
            //     projection of variable bounds - to be checked out once more 
            //----------------------------------------------------------------------------
            Vector lx =  (constraints_memory_.x_lower[level] - this->memory_.x[level]); 
            Scalar lower_multiplier = 1.0/constraints_memory_.P_inf_norm[level-1] * max(lx); 
            constraints_memory_.x_lower[level-1] = this->memory_.x[level-1] + local_values(local_size(this->memory_.x[level-1]).get(0), lower_multiplier);

            Vector ux =  (constraints_memory_.x_upper[level] - this->memory_.x[level]); 
            Scalar upper_multiplier = 1.0/constraints_memory_.P_inf_norm[level-1] * min(ux); 
            constraints_memory_.x_upper[level-1] = this->memory_.x[level-1] + local_values(local_size(this->memory_.x[level-1]).get(0), upper_multiplier);

            //----------------------------------------------------------------------------
            //     intersect bounds on coarse level
            //----------------------------------------------------------------------------

            constraints_memory_.active_upper[level-1] = local_zeros(local_size(this->memory_.x[level-1]).get(0)); 
            constraints_memory_.active_lower[level-1] = local_zeros(local_size(this->memory_.x[level-1]).get(0)); 
            {   
                Write<Vector>   rv(constraints_memory_.active_upper[level-1]), rw(constraints_memory_.active_lower[level-1]); 
                Read<Vector>    rl(this->memory_.x[level-1]), rq(constraints_memory_.x_lower[level-1]), re(constraints_memory_.x_upper[level-1]), rr(constraints_memory_.tr_lower[level-1]), rt(constraints_memory_.tr_upper[level-1]); 

                Range r = range(this->memory_.x[level-1]);

                for(SizeType i = r.begin(); i != r.end(); ++i)
                {
                    constraints_memory_.active_upper[level-1].set(i, std::min(constraints_memory_.tr_upper[level-1].get(i), constraints_memory_.x_upper[level-1].get(i))); 
                    constraints_memory_.active_lower[level-1].set(i, std::max(constraints_memory_.tr_lower[level-1].get(i), constraints_memory_.x_lower[level-1].get(i))); 
                }
            }
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
        virtual bool solve_qp_subproblem(const Matrix & H, const Vector & g, Vector & s, const SizeType & level, const bool & flg) override
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

            if(flg)
            {
                // setting should be really parameters from outside ... 
                this->_coarse_tr_subproblem->atol(1e-16); 
                this->_coarse_tr_subproblem->max_it(50); 

                if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->_coarse_tr_subproblem.get()))
                    tr_subproblem->tr_constrained_solve(H, g, s, box);
            }
            else
            {
                this->_smoother_tr_subproblem->atol(1e-16); 
                this->_smoother_tr_subproblem->max_it(1);
                
                if(TRSubproblem * tr_subproblem = dynamic_cast<TRSubproblem*>(this->_smoother_tr_subproblem.get()))
                    tr_subproblem->tr_constrained_solve(H, g, s, box);
            }

            return true; 
        }


    protected:   
        ConstraintsLevelMemory <Vector>         constraints_memory_;

        BoxConstraints box_constraints_;        // constraints on the finest level.... 
        bool has_box_constraints_;               // as we can run rmtr with inf. norm also without constraints... 


    };

}

#endif //UTOPIA_RMTR_HPP