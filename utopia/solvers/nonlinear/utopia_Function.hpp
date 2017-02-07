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

    };
}
#endif //UTOPIA_SOLVER_FUNCTION_HPP
