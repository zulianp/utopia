#ifndef UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP
#define UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{

    template<class Matrix, class Vector>
    class VariableBoundNonlinearSolver : public NonLinearSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;
        

        typedef utopia::LinearSolver<Matrix, Vector>            Solver;
        typedef utopia::BoxConstraints<Vector>                  BoxConstraints;
        
        
    public:
        VariableBoundNonlinearSolver( const std::shared_ptr <Solver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >(),
                                      const Parameters params = Parameters()):
        NonLinearSolver<Matrix, Vector>(linear_solver, params)
        {
            set_parameters(params);
        }
        
        bool solve(Function<Matrix, Vector> &fun, Vector &x) override = 0; 
       
       
        virtual void set_parameters(const Parameters params) override
        {
            NonLinearSolver<Matrix, Vector>::set_parameters(params);
        }
        
        
        virtual bool set_box_constraints(const BoxConstraints & box)
        {
            constraints_ = box;
            return true;
        }

        virtual bool  get_box_constraints(BoxConstraints & box)
        {
            box = constraints_;
            return true;
        }

        virtual const Vector & get_upper_bound()
        {
          if(constraints_.has_upper_bound())
            return *constraints_.upper_bound(); 
          else
          {
            utopia_error("upper bound does not exist. \n"); 
            return local_values(1, 9e9); 
          }

        }

        virtual const Vector & get_lower_bound()
        {
          if(constraints_.has_lower_bound())
            return *constraints_.lower_bound(); 
          else
          {
            utopia_error("lower bound does not exist. \n"); 
            return local_values(1, -9e9); 
          }
        }        


    protected: 
      virtual Scalar criticality_measure_infty(const Vector & x, const Vector & g)
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


      bool get_projection(const Vector & x, const Vector &lb, const Vector &ub, Vector & Pc)
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

    void make_iterate_feasible(Vector & x)
    {
        if(!constraints_.has_upper_bound() || !constraints_.has_lower_bound())
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


    // TO BE Changed...
    protected:
        BoxConstraints                  constraints_;

        
    };



}
#endif //UTOPIA_VARIABLE_BOUND_NONLINEAR_SOLVER_HPP
