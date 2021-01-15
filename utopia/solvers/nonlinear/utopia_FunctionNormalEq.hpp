#ifndef UTOPIA_SOLVER_NORMAL_EQ_FUNCTION_HPP
#define UTOPIA_SOLVER_NORMAL_EQ_FUNCTION_HPP

#include "utopia_Base.hpp"

namespace utopia {
    template <class Vector>
    class LeastSquaresFunctionBase {
    public:
        DEF_UTOPIA_SCALAR(Vector);

        virtual ~LeastSquaresFunctionBase() = default;

        // abstraction for normal equations
        virtual bool residual(const Vector & /*point*/, Vector & /*result*/) const = 0;
        virtual bool value(const Vector &, Scalar &val) const = 0;
        virtual bool update(const Vector &) { return true; }
    };

    /**
     * @brief      Nonlinear Function for normal equation problems.
     *             Additionally routines from Function class, one should supply also
     * residual and jacobian.
     *
     * @todo       Intorduce approximate Hessian updates strategies, e.g. BFGS, ...
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector>
    class LeastSquaresFunction : public LeastSquaresFunctionBase<Vector> {
    public:
        DEF_UTOPIA_SCALAR(Vector);

        ~LeastSquaresFunction() override = default;

        // abstraction for normal equations
        virtual bool jacobian(const Vector & /*x*/, Matrix & /*hessian*/) const = 0;
    };
}  // namespace utopia
#endif  // UTOPIA_SOLVER_NORMAL_EQ_FUNCTION_HPP