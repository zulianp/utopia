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
    class RMTR_inf :   public RMTR<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                    SizeType;

        // pay attention that this one is inf norm... 
        typedef utopia::TRBoxSubproblem<Matrix, Vector>     TRSubproblem; 

        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef utopia::Level<Matrix, Vector>               Level;


        typedef typename NonlinearMultiLevelBase<Matrix, Vector>::Fun Fun;


        typedef utopia::BoxConstraints<Vector>      BoxConstraints;
        typedef utopia::RMTR<Matrix, Vector>        RMTR;

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
            return "RMTR in infinity norm";  
        }
        

        /**
         * @brief      Sets the box constraints.
         *
         * @param      box   The box
         *
         */
        virtual bool set_box_constraints(BoxConstraints & box)
        {
          box_constraints_ = box; 
          has_box_constraints_ = true; 
          return true; 
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

            // on the finest level, no intersection with tr bounds is needed... 
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
        }



        virtual void init_coarse_level_constrains(const SizeType & level) override
        {
            std::cout<<"-------- to be done \n"; 
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
            std::cout<<"-------- to be done \n"; 

            return false; 
        }



        /**
         * @brief      "Heuristics", which decides if it makes sense to go to the coarse level or no
         *
         * @param[in]  g_restricted  Restricted gradient
         * @param[in]  g_coarse      Coarse level gradient 
         *
         */
        virtual bool grad_smoothess_termination(const Vector & g_restricted, const Vector & g_coarse) override
        {
            std::cout<<"-------- to be done \n";  
            return false;
        }


        virtual Scalar criticality_measure(const SizeType & level) override
        {
            std::cout<<"-------- to be done \n";  
            return 0.0;
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
            // this params should not be as hardcodded as they are...
            // if(flg)
            // {
            //     _coarse_tr_subproblem->atol(1e-16); 
            //     _coarse_tr_subproblem->max_it(5000); 
            //     _coarse_tr_subproblem->tr_constrained_solve(H, g, s, memory_.delta[level]); 
            // }
            // else
            // {
            //     _smoother_tr_subproblem->atol(1e-16); 
            //     _smoother_tr_subproblem->max_it(5);
            //     _smoother_tr_subproblem->tr_constrained_solve(H, g, s, memory_.delta[level]); 
            // }

            std::cout<<"----- to be done ------ \n"; 

            return true; 
        }



// .............................................. rmtr inf functions .............................................
        // correct only if P has only positive elements ... 
        // x is at level l and we are going to initialialize one level down... 
        void set_tr_constrain(const Vector &x, const Scalar level)
        {
            // Vector lower = local_zeros(local_size(x)), upper = local_zeros(local_size(x)); 
            // Vector upper_coarse, lower_coarse;
            
            // Vector tr_fine_last_lower = x - local_values(local_size(x).get(0), this->get_delta(level));   
            // Vector tr_fine_last_upper = x + local_values(local_size(x).get(0), this->get_delta(level));   

            // if(level < this->n_levels())
            // {              
            //     {   
            //         Read<Vector> rv(tr_fine_last_lower); 
            //         Read<Vector> rl(level_constraints_[level].tr_lower); 
            //         Write<Vector> wv(lower); 

            //         Range r = range(lower);

            //         for(SizeType i = r.begin(); i != r.end(); ++i) 
            //             lower.set(i, std::max(level_constraints_[level].tr_lower.get(i), tr_fine_last_lower.get(i))); 
            //     }
                
            //     {   
            //         Read<Vector> rv(tr_fine_last_upper); 
            //         Read<Vector> rl(level_constraints_[level].tr_upper); 
            //         Write<Vector> wv(upper); 

            //         Range r = range(upper);

            //         for(SizeType i = r.begin(); i != r.end(); ++i) 
            //             upper.set(i, std::min(level_constraints_[level].tr_upper.get(i), tr_fine_last_upper.get(i))); 
            //     }

            //     this->transfer(level-1).restrict(upper, upper_coarse);
            //     this->transfer(level-1).restrict(lower, lower_coarse);
            // }
            // else
            // {
            //     this->transfer(level-1).restrict(tr_fine_last_upper, upper_coarse);
            //     this->transfer(level-1).restrict(tr_fine_last_lower, lower_coarse);
            // }

            // level_constraints_[level-1].tr_upper = upper_coarse; 
            // level_constraints_[level-1].tr_lower = lower_coarse; 


            std::cout<<"--- to be done ----- \n"; 
        }



    protected:   
        ConstraintsLevelMemory <Vector>         constraints_memory_;

        BoxConstraints box_constraints_;        // constraints on the finest level.... 

        bool has_box_constraints_; // as we can run rmtr with inf. norm also without constraints... 


    };

}

#endif //UTOPIA_RMTR_HPP