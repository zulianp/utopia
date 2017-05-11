#ifndef UTOPIA_SOLVER_FUNCTION_HPP
#define UTOPIA_SOLVER_FUNCTION_HPP

#include "utopia_Base.hpp"

namespace utopia 
{
    /**
     * @brief      Base class for Nonlinear Function. All application context needed by solver is usually provided inside of this functions.
     *             In optimization settings, user needs to supply value(energy), gradient, hessian.
     *
     * @todo       Intorduce approximate Hessian updates strategies, e.g. BFGS, ... 
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector>
    class Function 
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        virtual ~Function() { }

        // TODO:: this needs to be changed ! 
        Function(const Vector & x_init = local_zeros(1), const Vector & rhs = local_zeros(1)) : 
                _x_init(x_init),
                _rhs(rhs)
        {

        }

        virtual bool value(const Vector &/*point*/, Scalar &/*value*/) const = 0;
        virtual bool gradient(const Vector &/*point*/, Vector &/*result*/) const = 0;
        virtual bool hessian(const Vector &x, Matrix &H) const = 0;

        virtual bool hessian(const Vector &/*point*/, Matrix &/*result*/, Matrix &/*preconditioner*/) const 
        {
            return false;
        }

        virtual bool has_preconditioner() const {
            return false;
        }

        virtual bool update(const Vector &/*point*/) { return true; };


        virtual bool set_rhs(const Vector & rhs)
        {
            _rhs = rhs; 
            return true; 
        }

        virtual bool reset_rhs()
        {
            _rhs = local_zeros(local_size(_rhs)); 
            return true; 
        }


        virtual bool get_rhs( Vector & rhs)
        {
            rhs = _rhs; 
            return true; 
        }

        virtual bool has_rhs() const
        {
            return (empty(_rhs))? false : true; 
        }



        virtual bool get_boundary_values(Vector & x)
        {   
            x = _x_init; 
            return true; 
        }


        virtual bool boundary_values_init(const Vector &x_in)
        {
            _x_init = x_in; 
            return true; 
        }




     protected:
        Vector _x_init; 
        Vector _rhs;


    };
}
#endif //UTOPIA_SOLVER_FUNCTION_HPP
