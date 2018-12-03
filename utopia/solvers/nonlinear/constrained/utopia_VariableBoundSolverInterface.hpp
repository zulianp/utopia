#ifndef UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP
#define UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_BoxConstraints.hpp"


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
        VariableBoundSolverInterface()
        {

        }

        virtual ~VariableBoundSolverInterface()
        {

        }
                
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

        Vector Pc; 
        Vector x_g = x - g; 
        Vector ub, lb; 

        Scalar n = local_size(x).get(0); 

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


      bool get_projection(const Vector & x, const Vector &lb, const Vector &ub, Vector & Pc) const 
      {
          Pc = local_values(local_size(x).get(0), 1.0);
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

          return true;
      }

    void make_iterate_feasible(Vector & x) const 
    {
        if(!constraints_.has_upper_bound() && !constraints_.has_lower_bound())
            return; 

        const Vector x_old = x; 

        if(constraints_.has_upper_bound() && constraints_.has_lower_bound())
        {
            const auto &ub = *constraints_.upper_bound();
            const auto &lb = *constraints_.lower_bound();

            {
              Read<Vector> r_ub(ub), r_lb(lb), r_x(x_old);
              Write<Vector> wv(x); 

              each_write(x, [ub, lb, x_old](const SizeType i) -> double { 
                          Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x_old.get(i);  
                          if(li >= xi)
                            return li; 
                          else
                            return (ui <= xi) ? ui : xi; }   );
            }
        }
        else if(constraints_.has_upper_bound() && !constraints_.has_lower_bound())
        {
            const auto &ub = *constraints_.upper_bound();

            {
              Read<Vector> r_ub(ub), r_x(x_old);
              Write<Vector> wv(x); 

              each_write(x, [ub, x_old](const SizeType i) -> double { 
                          Scalar ui =  ub.get(i); Scalar xi =  x_old.get(i);  
                            return (ui <= xi) ? ui : xi; }   );
            }
        }
        else
        {
            const auto &lb = *constraints_.lower_bound();

            {
              Read<Vector> r_lb(lb), r_x(x_old);
              Write<Vector> wv(x); 

              each_write(x, [lb, x_old](const SizeType i) -> double { 
                          Scalar li =  lb.get(i); Scalar xi =  x_old.get(i);  
                          return (li >= xi) ? li : xi; }   );
            }
        }    

    }      


      virtual BoxConstraints  merge_pointwise_constraints_with_uniform_bounds(const Vector & x_k, const Scalar & lb_uniform, const Scalar & ub_uniform) const 
      {
          Vector l_f, u_f; 

          if(constraints_.has_upper_bound())
          {
              Vector u =  *constraints_.upper_bound() - x_k; 
              u_f = local_zeros(local_size(x_k).get(0)); 
              {   
                  Read<Vector> rv(u); 
                  Write<Vector> wv(u_f); 

                  each_write(u_f, [ub_uniform, u](const SizeType i) -> double { 
                      return  (u.get(i) <= ub_uniform)  ? u.get(i) : ub_uniform; }   );
              }
          }
          else
              u_f = local_values(local_size(x_k).get(0), ub_uniform); ; 

          if(constraints_.has_lower_bound())
          {
              Vector l = *(constraints_.lower_bound()) - x_k; 
              l_f = local_zeros(local_size(x_k).get(0)); 

              {   
                  Read<Vector> rv(l); 
                  Write<Vector> wv(l_f); 

                  each_write(l_f, [lb_uniform, l](const SizeType i) -> double { 
                      return  (l.get(i) >= lb_uniform)  ? l.get(i) : lb_uniform;  }   );
              }
          }
          else
              l_f = local_values(local_size(x_k).get(0), lb_uniform); 

          return make_box_constaints(std::make_shared<Vector>(l_f), std::make_shared<Vector>(u_f));
      }



    protected:
        BoxConstraints                  constraints_;

        
    };



}
#endif //UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP
