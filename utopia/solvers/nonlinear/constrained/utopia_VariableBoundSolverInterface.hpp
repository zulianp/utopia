#ifndef UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP
#define UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_BoxConstraints.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_For.hpp"


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
      virtual Scalar criticality_measure_infty(const Vector & x, const Vector & g) const
      {

        //FIXME remove temporaries
        Vector Pc;
        Vector x_g = x - g;
        Vector ub, lb;

        Scalar n = local_size(x);

        if(constraints_.has_upper_bound())
          ub = *constraints_.upper_bound();
        else
          ub =  local_values(n, 9e15);

        if(constraints_.has_lower_bound())
          lb = *constraints_.lower_bound();
        else
          lb =  local_values(n, -9e15);

        get_projection(x_g, lb, ub, Pc);

        Pc -= x;
        return norm2(Pc);
      }

  public:  // expose it for CUDA
      bool get_projection(const Vector & x, const Vector &lb, const Vector &ub, Vector & Pc) const
      {
        if(empty(Pc)){
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

        //FIXME remove temporaries
        const Vector x_old = x;

        if(constraints_.has_upper_bound() && constraints_.has_lower_bound())
        {
            const auto &ub = *constraints_.upper_bound();
            const auto &lb = *constraints_.lower_bound();

          {
            auto d_lb     = const_device_view(lb);
            auto d_ub     = const_device_view(ub);
            auto d_xold   = const_device_view(x_old);

            parallel_each_write(x, UTOPIA_LAMBDA(const SizeType i) -> Scalar
            {
                Scalar li = d_lb.get(i);
                Scalar ui = d_ub.get(i);
                Scalar xi = d_xold.get(i);

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
              auto d_xold   = const_device_view(x_old);

              parallel_each_write(x, UTOPIA_LAMBDA(const SizeType i) -> Scalar
              {
                Scalar ui = d_ub.get(i);
                Scalar xi = d_xold.get(i);
                return (ui <= xi) ? ui : xi;
              });
            }
        }
        else
        {
            const auto &lb = *constraints_.lower_bound();

            {
              auto d_lb     = const_device_view(lb);
              auto d_xold   = const_device_view(x_old);

              parallel_each_write(x, UTOPIA_LAMBDA(const SizeType i) -> Scalar
              {
                Scalar li =  d_lb.get(i);
                Scalar xi =  d_xold.get(i);
                return (li >= xi) ? li : xi;
              });
            }
        }

    }

      virtual BoxConstraints  merge_pointwise_constraints_with_uniform_bounds(const Vector & x_k, const Scalar & lb_uniform, const Scalar & ub_uniform) const
      {
          //FIXME remove temporaries
          Vector l_f, u_f;

          if(constraints_.has_upper_bound())
          {
              Vector u =  *constraints_.upper_bound() - x_k;
              u_f = local_zeros(local_size(x_k));

              {
                auto d_u    = const_device_view(u);

                parallel_each_write(u_f, UTOPIA_LAMBDA(const SizeType i) -> Scalar
                {
                  auto val = d_u.get(i);
                  return  (val <= ub_uniform)  ? val : ub_uniform;
                });
              }
          }
          else
              u_f = local_values(local_size(x_k), ub_uniform); ;

          if(constraints_.has_lower_bound())
          {
              Vector l = *(constraints_.lower_bound()) - x_k;
              l_f = local_zeros(local_size(x_k));

            {
              auto d_l    = const_device_view(l);

              parallel_each_write(l_f, UTOPIA_LAMBDA(const SizeType i) -> Scalar
              {
                  auto val = d_l.get(i);
                  return  (val >= lb_uniform)  ? val : lb_uniform;
              });
            }



          }
          else
              l_f = local_values(local_size(x_k), lb_uniform);

          return make_box_constaints(std::make_shared<Vector>(l_f), std::make_shared<Vector>(u_f));
      }


      virtual BoxConstraints  build_correction_constraints(const Vector & x_k) const
      {
          //FIXME remove temporaries
          Vector l_f, u_f;

          if(constraints_.has_upper_bound())
              u_f =  *constraints_.upper_bound() - x_k;
          else
              u_f = local_values(local_size(x_k), 9e12);

          if(constraints_.has_lower_bound())
              l_f = *(constraints_.lower_bound()) - x_k;
          else
              l_f = local_values(local_size(x_k), -9e12);

          return make_box_constaints(std::make_shared<Vector>(l_f), std::make_shared<Vector>(u_f));
      }


    protected:
        BoxConstraints                  constraints_;


    };

}
#endif //UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP
