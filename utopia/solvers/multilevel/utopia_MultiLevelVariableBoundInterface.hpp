#ifndef UTOPIA_MULTILEVEL_VARIABLE_BOUND_INTERFACE_HPP
#define UTOPIA_MULTILEVEL_VARIABLE_BOUND_INTERFACE_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_LevelMemory.hpp"

#include "utopia_IdentityTransfer.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{
    template<class Vector>
    class MultilevelVariableBoundSolverInterface 
    {
        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

        typedef utopia::BoxConstraints<Vector>                  BoxConstraints;


    public:
        MultilevelVariableBoundSolverInterface(): has_box_constraints_(false)
        {

        }

        virtual ~MultilevelVariableBoundSolverInterface() = default;


        virtual void set_box_constraints(BoxConstraints & box)
        {
          box_constraints_ = box;
          has_box_constraints_ = true;
        }

        virtual BoxConstraints & get_box_constraints()
        {
          return box_constraints_;
        }


      protected:
        void init_memory(const std::vector<SizeType> & n_dofs_)
        {
            constraints_memory_.init_memory(n_dofs_);

            SizeType fine_local_size = n_dofs_.back(); 
            SizeType n_levels =  n_dofs_.size(); 

            const SizeType fine_level =n_levels -1;
            const Scalar inf = std::numeric_limits<Scalar>::infinity();

            if(has_box_constraints_)
            {
                if(box_constraints_.has_upper_bound()){
                    constraints_memory_.x_upper[fine_level] = *box_constraints_.upper_bound();
                }
                else{
                    constraints_memory_.x_upper[fine_level].set(inf); 
                }

                if(box_constraints_.has_lower_bound()){
                    constraints_memory_.x_lower[fine_level] = *box_constraints_.lower_bound();
                }
                else{
                    constraints_memory_.x_lower[fine_level].set(-1.0 * inf);
                }

                constraints_memory_.active_upper[fine_level] = constraints_memory_.x_upper[fine_level];
                constraints_memory_.active_lower[fine_level] = constraints_memory_.x_lower[fine_level];
            }
            else
            {
                constraints_memory_.active_upper[fine_level].set(inf);
                constraints_memory_.active_lower[fine_level].set(-1.0 * inf);                
            }

            // inherited tr bound constraints...
            constraints_memory_.tr_upper[fine_level].set(inf);
            constraints_memory_.tr_lower[fine_level].set(-1.0 * inf);            


            // // precompute norms of prolongation operators needed for projections of constraints...
            // for(auto l = 0; l < fine_level; l++)
            //     constraints_memory_.P_inf_norm[l] = this->transfer(l).interpolation_inf_norm();

        }



        virtual bool check_feasibility(const SizeType & level, const Vector & x)
        {
            bool terminate = false;

            {
                Read<Vector> ru(constraints_memory_.tr_upper[level]);
                Read<Vector> rl(constraints_memory_.tr_lower[level]);
                Read<Vector> rx(x);

                Range r = range(constraints_memory_.tr_upper[level]);

                for(SizeType i = r.begin(); i != r.end(); ++i)
                {
                    Scalar xi = x.get(i);
                    Scalar li = constraints_memory_.tr_lower[level].get(i);
                    Scalar ui = constraints_memory_.tr_upper[level].get(i);

                   if(xi < li || xi > ui)
                        terminate = true;
                }
            }


           //  using ForLoop = utopia::ParallelFor<Traits<Vector>::Backend>;

           // {
           //      auto d_u = const_device_view(constraints_memory_.tr_upper[level]);
           //      auto d_l = const_device_view(constraints_memory_.tr_lower[level]);
           //      auto d_x = const_device_view(x);

           //      ForLoop::apply(range(x), UTOPIA_LAMBDA(const SizeType i)
           //      {
           //          const Scalar xi = d_x.get(i);
           //          const Scalar li = d_l.get(i);
           //          const Scalar ui = d_u.get(i);

           //          // if(xi < li || xi > ui){
           //          //     terminate = true;
           //          // }

           //      });
           //  }

            return x.comm().disjunction(terminate);
          
        }


      void get_projection(const Vector &lb, const Vector &ub, Vector & x)
      {

        {
            auto d_lb = const_device_view(lb);
            auto d_ub = const_device_view(ub);

            parallel_transform(x, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar 
            {
                Scalar li = d_lb.get(i);
                Scalar ui = d_ub.get(i);
                if(li >= xi){
                  return li;
                }
                else{
                  return (ui <= xi) ? ui : xi;
                }
            });
        }

      }        


        virtual Scalar criticality_measure_inf(const SizeType & level, const Vector & x, const Vector & g)
        {
            constraints_memory_.help[level] = x - g; 
            get_projection(constraints_memory_.active_lower[level], constraints_memory_.active_upper[level], constraints_memory_.help[level]);
            constraints_memory_.help[level] -= x;

            return norm2(constraints_memory_.help[level]);
        }




        template<typename Matrix>
        void init_level(const SizeType & level, const Vector & x_finer_level,  const Vector & x_level, const Scalar & delta_fine, Transfer<Matrix, Vector> & transfer)
        {
            const SizeType finer_level = level+1;

            if( IdentityTransfer<Matrix, Vector>* id_transfer =  dynamic_cast<IdentityTransfer<Matrix, Vector>* > (&transfer))
            {
                constraints_memory_.tr_upper[level] =  x_finer_level + local_values(local_size(x_finer_level), delta_fine);
                constraints_memory_.tr_lower[level] =  x_finer_level - local_values(local_size(x_finer_level), delta_fine);
            }
            else
            {
                //----------------------------------------------------------------------------
                //     soft projection of tr bounds
                //----------------------------------------------------------------------------
                Vector tr_fine_last_lower = x_finer_level - local_values(local_size(x_finer_level), delta_fine);
                {
                    ReadAndWrite<Vector> rv(tr_fine_last_lower);
                    Read<Vector> rl(constraints_memory_.tr_lower[finer_level]);

                    Range r = range(tr_fine_last_lower);

                    for(SizeType i = r.begin(); i != r.end(); ++i){
                        tr_fine_last_lower.set(i, std::max(constraints_memory_.tr_lower[finer_level].get(i), tr_fine_last_lower.get(i)));
                    }
                }

                Vector tr_fine_last_upper = x_finer_level + local_values(local_size(x_finer_level), delta_fine);
                {
                    ReadAndWrite<Vector> rv(tr_fine_last_upper);
                    Read<Vector> rl(constraints_memory_.tr_upper[finer_level]);

                    Range r = range(tr_fine_last_upper);

                    for(SizeType i = r.begin(); i != r.end(); ++i){
                        tr_fine_last_upper.set(i, std::min(constraints_memory_.tr_upper[finer_level].get(i), tr_fine_last_upper.get(i)));
                    }
                }

                //------------------------ new version, taking into account  positive and negative elements projection separatelly -----------------
                transfer.project_down_positive_negative(tr_fine_last_lower, tr_fine_last_upper, constraints_memory_.tr_lower[level]);
                transfer.project_down_positive_negative(tr_fine_last_upper, tr_fine_last_lower, constraints_memory_.tr_upper[level]);

            }


            if(has_box_constraints_)
            {

                if( IdentityTransfer<Matrix, Vector>* id_transfer =  dynamic_cast<IdentityTransfer<Matrix, Vector>* > (&transfer))
                {
                    constraints_memory_.x_lower[level] = constraints_memory_.x_lower[finer_level]; 
                    constraints_memory_.x_upper[level] = constraints_memory_.x_upper[finer_level]; 
                }
                else
                {
                    //----------------------------------------------------------------------------
                    //     projection of variable bounds to the coarse level
                    //----------------------------------------------------------------------------
                    Vector lx =  (constraints_memory_.x_lower[finer_level] - x_finer_level);
                    Scalar lower_multiplier = 1.0/constraints_memory_.P_inf_norm[level] * max(lx);
                    constraints_memory_.x_lower[level] = x_level + local_values(local_size(x_level), lower_multiplier);

                    Vector ux =  (constraints_memory_.x_upper[finer_level] - x_finer_level);
                    Scalar upper_multiplier = 1.0/constraints_memory_.P_inf_norm[level] * min(ux);
                    constraints_memory_.x_upper[level] = x_level + local_values(local_size(x_level), upper_multiplier);
                }

                //----------------------------------------------------------------------------
                //     intersect bounds on the coarse level
                //----------------------------------------------------------------------------
                constraints_memory_.active_upper[level] = local_zeros(local_size(x_level));
                constraints_memory_.active_lower[level] = local_zeros(local_size(x_level));
                {
                    Write<Vector>   rv(constraints_memory_.active_upper[level]), rw(constraints_memory_.active_lower[level]);
                    Read<Vector>    rl(x_level), rq(constraints_memory_.x_lower[level]), re(constraints_memory_.x_upper[level]), rr(constraints_memory_.tr_lower[level]), rt(constraints_memory_.tr_upper[level]);

                    Range r = range(x_level);

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
        }



    protected:
        BoxConstraints box_constraints_;        // constraints on the finest level....
        bool has_box_constraints_;               // as we can run rmtr with inf. norm also without constraints...      

        ConstraintsLevelMemory <Vector>         constraints_memory_;


    };

}
#endif //UTOPIA_MULTILEVEL_VARIABLE_BOUND_INTERFACE_HPP
