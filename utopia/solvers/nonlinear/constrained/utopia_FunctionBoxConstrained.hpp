#ifndef UTOPIA_FUNCTION_BOX_CONSTRAINTS_HPP
#define UTOPIA_FUNCTION_BOX_CONSTRAINTS_HPP

#include "utopia_Base.hpp"

namespace utopia
{
    /**
     * @brief       Function for box-constrained nonlinear problems.
     *              Additionally to gradient, hessian, function value provides interface to lower and upper bound on solution.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector>
    class FunctionBoxConstrained : public Function<Matrix, Vector>
    {
        public:
            DEF_UTOPIA_SCALAR(Matrix)

            virtual ~FunctionBoxConstrained() { }

            virtual bool upper_bound(Vector &/*ub*/) const = 0;
            virtual bool lower_bound(Vector &/*lb*/) const = 0;

    };
}
#endif //UTOPIA_FUNCTION_BOX_CONSTRAINTS_HPP
