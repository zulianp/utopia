#ifndef UTOPIA_EXTENDED_FUNCTION_HPP
#define UTOPIA_EXTENDED_FUNCTION_HPP

#include "utopia_Base.hpp"

namespace utopia 
{
    /**
     * @brief      Class for Nonlinear Function, all application context needed by solver is usually provided inside of this functions.
     *             In optimization settings, user needs to supply value(energy), gradient, hessian.
     *             Difference, between Function and ExtendedFunction is that here, we make use of additional informations to improve convergence
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector>
    class ExtendedFunction : public Function<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        virtual ~ExtendedFunction() { }


        ExtendedFunction(const Vector & x_init, const Vector & bc_marker, const Vector & rhs) : 
                _x_eq_values(x_init),
                _rhs(rhs), 
                _eq_constrains_flg(bc_marker)
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

        virtual bool get_eq_constrains_values(Vector & x)
        {   
            x = _x_eq_values; 
            return true; 
        }

        virtual bool get_eq_constrains_flg(Vector & x)
        {   
            x = _eq_constrains_flg; 
            return true; 
        }

        virtual bool set_equality_constrains(const Vector &eq_constrains_flg, const Vector &x_in)
        {
            _x_eq_values             =  x_in; 
            _eq_constrains_flg  = eq_constrains_flg; 
            return true; 
        }

     protected:
        Vector _x_eq_values; 
        Vector _rhs;
        Vector _eq_constrains_flg;


    };
}
#endif //UTOPIA_EXTENDED_FUNCTION_HPP
