#ifndef UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP
#define UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_For.hpp"
#include "utopia_Allocations.hpp"


#include <iomanip>
#include <limits>


namespace utopia
{
    template<class Vector>
    class VariableBoundSolverInterface
    {
        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

        typedef utopia::BoxConstraints<Vector>                  BoxConstraints;


    public:
        VariableBoundSolverInterface() = default;
        virtual ~VariableBoundSolverInterface() = default;

        virtual bool set_box_constraints(const BoxConstraints & box)
        {
            constraints_ = box;
            return true;
        }

        virtual const BoxConstraints & get_box_constraints() const
        {
            return constraints_;
        }

        virtual const Vector & get_upper_bound() const
        {
          if(!constraints_.has_upper_bound()){
            utopia_error("VariableBoundSolverInterface::upper bound does not exist. \n");
          }

          return *constraints_.upper_bound();
        }

        virtual const Vector & get_lower_bound() const
        {
          if(!constraints_.has_lower_bound()){
            utopia_error("VariableBoundSolverInterface::lower bound does not exist. \n");
          }

          return *constraints_.lower_bound();
        }

        virtual bool has_bound() const
        {
          return constraints_.has_bound();
        }

        virtual bool has_lower_bound() const
        {
          return constraints_.has_lower_bound();
        }

        virtual bool has_upper_bound() const
        {
          return constraints_.has_upper_bound();
        }

        virtual void fill_empty_bounds()
        {
          return constraints_.fill_empty_bounds();
        }



    protected:
      virtual Scalar criticality_measure_infty(const Vector & x, const Vector & g)
      {
        if(empty(Pc_) || size(Pc_) != size(x)){
          Pc_ = 0.0 * x; 
        }

        xg_ = x - g;

        if(!constraints_.has_upper_bound() || !constraints_.has_lower_bound())
        {
          this->fill_empty_bounds(); 
        }

        const auto &ub = *constraints_.upper_bound();
        const auto &lb = *constraints_.lower_bound();

        get_projection(xg_, lb, ub, Pc_);
        Pc_ -= x;

        return norm2(Pc_);
      }

  public:  // expose it for CUDA
      bool get_projection(const Vector & x, const Vector &lb, const Vector &ub, Vector & Pc) const
      {
        if(empty(Pc) || size(Pc) != size(x)){
            Pc = 0.0*x;
        }

        {
            auto d_lb = const_device_view(lb);
            auto d_ub = const_device_view(ub);
            auto d_x  = const_device_view(x);

            parallel_each_write(Pc, UTOPIA_LAMBDA(const SizeType i) -> Scalar
            {
                Scalar li = d_lb.get(i);
                Scalar ui = d_ub.get(i);
                Scalar xi = d_x.get(i);

                if(li >= xi)
                  return li;
                else
                  return (ui <= xi) ? ui : xi;

            });
        }

        return true;
      }

    void make_iterate_feasible(Vector & x) const
    {
        if(!constraints_.has_upper_bound() && !constraints_.has_lower_bound())
            return;

        if(constraints_.has_upper_bound() && constraints_.has_lower_bound())
        {
            const auto &ub = *constraints_.upper_bound();
            const auto &lb = *constraints_.lower_bound();

          {
            auto d_lb     = const_device_view(lb);
            auto d_ub     = const_device_view(ub);

          parallel_transform(
                          x,
                          UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar {
                            Scalar li = d_lb.get(i);
                            Scalar ui = d_ub.get(i);
                            if(li >= xi)
                              return li;
                            else
                              return (ui <= xi) ? ui : xi;
                      });

          }

        }
        else if(constraints_.has_upper_bound() && !constraints_.has_lower_bound())
        {
            const auto &ub = *constraints_.upper_bound();

            {
              auto d_ub     = const_device_view(ub);

              parallel_transform(
                              x,
                              UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar {
                                Scalar ui = d_ub.get(i);
                                return (ui <= xi) ? ui : xi;
                          });
            }
        }
        else
        {
            const auto &lb = *constraints_.lower_bound();

            {
              auto d_lb     = const_device_view(lb);

              parallel_transform(
                              x,
                              UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar {
                                Scalar li =  d_lb.get(i);
                                return (li >= xi) ? li : xi;
                          });

            }
        }

    }

      virtual const BoxConstraints &  merge_pointwise_constraints_with_uniform_bounds(const Vector & x_k, const Scalar & lb_uniform, const Scalar & ub_uniform) 
      {
        correction_constraints_.fill_empty_bounds(local_size(x_k)); 

        if(constraints_.has_upper_bound())
        {
            auto &ub_merged = *correction_constraints_.upper_bound();
            ub_merged =  *constraints_.upper_bound() - x_k;

          {
            parallel_transform(ub_merged,
                              UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar 
                              {
                                return  (xi <= ub_uniform)  ? xi : ub_uniform; 
                              });
          }

        }
        else{
          correction_constraints_.upper_bound()->set(ub_uniform); 
        }


        if(constraints_.has_lower_bound())
        {
            auto &lb_merged = *correction_constraints_.lower_bound();
            lb_merged =  *constraints_.lower_bound() - x_k;

          {
            parallel_transform(lb_merged,
                              UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar 
                              {
                                return  (xi >= lb_uniform)  ? xi : lb_uniform; 
                              });
          }

        }
        else{
          correction_constraints_.lower_bound()->set(lb_uniform); 
        }

        return correction_constraints_; 
      }


      virtual const BoxConstraints & build_correction_constraints(const Vector & x_k) 
      {
          correction_constraints_.fill_empty_bounds(local_size(x_k)); 

          if(constraints_.has_upper_bound()){
            *correction_constraints_.upper_bound() =  *constraints_.upper_bound() - x_k;
          }
          else{
            correction_constraints_.upper_bound()->set(9e12); 
          }


          if(constraints_.has_lower_bound()){
            *correction_constraints_.lower_bound() = *(constraints_.lower_bound()) - x_k;
          }
          else{
            correction_constraints_.lower_bound()->set(-9e12); 
          }

          return correction_constraints_; 
      }


    protected:
        BoxConstraints                  constraints_;             // variable bound constraints 
        BoxConstraints                  correction_constraints_;  // constraints needed for correction 

        Vector Pc_;
        Vector xg_; 


    };

}
#endif //UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP
