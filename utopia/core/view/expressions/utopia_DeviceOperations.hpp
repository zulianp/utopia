#ifndef UTOPIA_DEVICE_OPERATIONS_HPP
#define UTOPIA_DEVICE_OPERATIONS_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_Operators.hpp"
#include "utopia_DeviceNorm.hpp"

#ifdef WITH_TRILINOS
#include <Kokkos_ArithTraits.hpp>
#endif

#include <limits>
#include <algorithm>
#include <cmath>

namespace utopia {

    template<class Expr>
    UTOPIA_INLINE_FUNCTION DeviceBinary<Expr, DeviceNumber<double>, Minus> operator-(const DeviceExpression<Expr> &left, const double &right)
    {
        return DeviceBinary<Expr, DeviceNumber<double>, Minus>(left, right);
    }

    template<class Expr>
    UTOPIA_INLINE_FUNCTION double operator-(const DeviceNumber<double> &left, const DeviceNumber<double> &right)
    {
        return static_cast<double>(left) - static_cast<double>(right);
    }

    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceNegate<Derived> operator-(const DeviceExpression<Derived> &expr) {
        return DeviceNegate<Derived>(expr.derived());
    }

    template<class Left, class Right>
    UTOPIA_INLINE_FUNCTION DeviceBinary<Left, Right, Minus> operator-(const DeviceExpression<Left> &left, const DeviceExpression<Right> &right) {
        return DeviceBinary<Left, Right, Minus>(left.derived(), right.derived());
    }

    template<class Left, class Right>
    UTOPIA_INLINE_FUNCTION DeviceBinary<Left, Right, Plus> operator+(const DeviceExpression<Left> &left, const DeviceExpression<Right> &right) {
        return DeviceBinary<Left, Right, Plus>(left.derived(), right.derived());
    }

    //Switch left with right
    template<class Left, class L, class R>
    UTOPIA_INLINE_FUNCTION DeviceBinary<DeviceMultiply<TensorView<L, 2>, TensorView<R, 1>>, TensorView<Left, 1>, Plus> operator+(
        const TensorView<Left, 1> &left,
        const DeviceMultiply<TensorView<L, 2>, TensorView<R, 1>> &right) {
        return DeviceBinary<DeviceMultiply<TensorView<L, 2>, TensorView<R, 1>>, TensorView<Left, 1>, Plus>(right.derived(), left.derived());
    }

    template<class Left, class Right>
    UTOPIA_INLINE_FUNCTION DeviceMultiply<Left, Right> operator*(
        const DeviceExpression<Left> &left,
        const DeviceExpression<Right> &right
    ) {
        return DeviceMultiply<Left, Right>(left.derived(), right.derived());
    }

    template<class Left>
    UTOPIA_INLINE_FUNCTION DeviceBinary<
        DeviceNumber<typename Left::Scalar>,
        Left,
        Multiplies
    >
    operator*(
        const DeviceExpression<Left> &left,
        const typename Left::Scalar &right
    ) {
        return DeviceBinary<DeviceNumber<typename Left::Scalar>, Left, Multiplies>(right, left.derived());
    }

    template<class Right>
    UTOPIA_INLINE_FUNCTION DeviceBinary<
        DeviceNumber<typename Right::Scalar>,
        Right,
        Multiplies
    >
    operator*(
        const typename Right::Scalar &left,
        const DeviceExpression<Right> &right
    ) {
        return DeviceBinary<DeviceNumber<typename Right::Scalar>, Right, Multiplies>(left, right.derived());
    }

    template<class Derived, class Operation>
    UTOPIA_INLINE_FUNCTION DeviceUnary<Derived, Operation> transform(const DeviceExpression<Derived> &expr, const Operation operation = Operation()) {
        return DeviceUnary<Derived, Operation>(expr, operation);
    }

    /**     @defgroup   transforms Transforms
     *       @ingroup    algebra
    */


    /**
     * @ingroup transforms
     * @brief   \f$ x_i^2  \f$.
     */
    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceUnary<Derived, Pow2> pow2(const DeviceExpression<Derived> &expr)
    {
        return DeviceUnary<Derived, Pow2>(expr.derived());
    }

    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceUnary<Derived, Log> logn(const DeviceExpression<Derived> &expr)
    {
        return DeviceUnary<Derived, Log>(expr.derived());
    }

    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceUnary<Derived, Exp> exp(const DeviceExpression<Derived> &expr)
    {
        return DeviceUnary<Derived, Exp>(expr.derived());
    }

    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceUnary<Derived, Cos> cos(const DeviceExpression<Derived> &expr)
    {
        return DeviceUnary<Derived, Cos>(expr.derived());
    }

    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceUnary<Derived, Sin> sin(const DeviceExpression<Derived> &expr)
    {
        return DeviceUnary<Derived, Sin>(expr.derived());
    }



    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceUnary<Derived, Pow> power(
        const DeviceExpression<Derived> &expr,
        const typename Traits<Derived>::Scalar &a)
    {
        return DeviceBinary<
            Derived,
            DeviceNumber<typename Traits<Derived>::Scalar>,
            Pow
            >(expr.derived(), a);
    }

    /**
     * @ingroup transforms
     * @brief   \f$ | x_i |  \f$.
     */
    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceUnary<Derived, Abs> abs(const DeviceExpression<Derived> &expr)
    {
        return DeviceUnary<Derived, Abs>(expr.derived());
    }

    /**
     * @ingroup transforms
     * @brief   \f$ \sqrt(x_i)  \f$.
     */
    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceUnary<Derived, Sqrt> sqrt(const DeviceExpression<Derived> &expr)
    {
        return DeviceUnary<Derived, Sqrt>(expr.derived());
    }

    /**
     * @ingroup transforms
     * @brief   \f$ \frac 1 {x_i}  \f$.
     */
    template<class Derived>
    UTOPIA_INLINE_FUNCTION
    DeviceBinary<
        DeviceNumber<typename Derived::Scalar>,
        Derived,
        Divides
    > operator / (const typename Derived::Scalar &left,
                  const DeviceExpression<Derived> &expr)
    {
        return DeviceBinary<
            DeviceNumber<typename Derived::Scalar>,
            Derived,
            Divides
        >(left, expr.derived());
    }

    template<class LDerived, class RDerived>
    UTOPIA_INLINE_FUNCTION DeviceBinary<LDerived, RDerived, Divides> operator / (
        const DeviceExpression<LDerived> &left,
        const DeviceExpression<RDerived> &right
        )
    {
        return DeviceBinary<LDerived, RDerived, Divides>(left.derived(), right.derived());
    }

    template<class Left, class Right>
    UTOPIA_INLINE_FUNCTION bool approxeq(const DeviceExpression<Left> &left,
                                         const DeviceExpression<Right> &right,
                                         const typename Right::Scalar tol = 1e-6)
    {
        return DeviceApproxEqual<Left, Right>::apply(left.derived(), right.derived(), tol);
    }

    template<class Left, class Right>
    UTOPIA_INLINE_FUNCTION typename Traits<Right>::Scalar dot(
        const DeviceExpression<Left> &left,
        const DeviceExpression<Right> &right)
    {
        return DeviceDot<Left, Right, Traits<Left>::Order>::apply(left.derived(), right.derived());
    }

    template<class Left, class Right>
    UTOPIA_INLINE_FUNCTION typename Traits<Right>::Scalar inner(
        const DeviceExpression<Left> &left,
        const DeviceExpression<Right> &right)
    {
        return dot(left.derived(), right.derived());
    }

    template<class Left, class Right>
    UTOPIA_INLINE_FUNCTION typename Traits<Right>::Scalar inner(
        const DeviceExpression<Left> &left,
        const typename Traits<Right>::Scalar &right)
    {
        return sum(left * right);
    }


    template<class Left, class Right>
    UTOPIA_INLINE_FUNCTION DeviceOuterProduct<Left, Right> outer(
        const DeviceExpression<Left> &left,
        const DeviceExpression<Right> &right)
    {
        return DeviceOuterProduct<Left, Right>(left.derived(), right.derived());
    }

    UTOPIA_INLINE_FUNCTION constexpr double inner(const double &left, const double &right)
    {
        return left * right;
    }

    UTOPIA_INLINE_FUNCTION constexpr float inner(const float &left, const float &right)
    {
        return left * right;
    }

    template<class Expr>
    UTOPIA_INLINE_FUNCTION typename Traits<Expr>::Scalar trace(
        const DeviceExpression<Expr> &expr)
    {
        return DeviceTrace<Expr>::apply(expr.derived());
    }

    template<class Expr>
    UTOPIA_INLINE_FUNCTION typename Traits<Expr>::Scalar norm1(
        const DeviceExpression<Expr> &expr)
    {
        return DeviceNorm<Expr, 1>::apply(expr.derived());
    }

    template<class Expr>
    UTOPIA_INLINE_FUNCTION typename Traits<Expr>::Scalar norm2(
        const DeviceExpression<Expr> &expr)
    {
        return DeviceNorm<Expr, 2>::apply(expr.derived());
    }

    template<class Expr>
    UTOPIA_INLINE_FUNCTION typename Traits<Expr>::Scalar norm_infty(
        const DeviceExpression<Expr> &expr)
    {
        return DeviceNorm<Expr, INFINITY_NORM_TAG>::apply(expr.derived());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise multiplication.
     */
    template<class Left, class Right>
    UTOPIA_INLINE_FUNCTION DeviceBinary<Left, Right, EMultiplies> e_mul(const DeviceExpression<Left> &left, const DeviceExpression<Right> &right) {
        return DeviceBinary<Left, Right, EMultiplies>(left.derived(),
                                                right.derived());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise min.
     */
    template<class Left, class Right>
    UTOPIA_INLINE_FUNCTION DeviceBinary<Left, Right, Min> min(const DeviceExpression<Left> &left, const DeviceExpression<Right> &right) {
        return DeviceBinary<Left, Right, Min>(left.derived(), right.derived());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise min.
     */
    // template<class Left, class Right, int Order>
    // UTOPIA_INLINE_FUNCTION DeviceBinary<Left, DeviceNumber<Right>, Min> min(const DeviceExpression<Left> &left, const Factory<Values<Right>, Order> &right) {
    //     return DeviceBinary<Left, DeviceNumber<Right>, Min>(left.derived(), right.type().value());
    // }

    /**
     * @ingroup elementwise
     * @brief   Pointwise min.
     */
    // template<class Left, class Right, int Order>
    // UTOPIA_INLINE_FUNCTION DeviceBinary<Left, DeviceNumber<Right>, Min> min(const Factory<Values<Right>, Order> &right, const DeviceExpression<Left> &left) {
    //     return DeviceBinary<Left, DeviceNumber<Right>, Min>(left.derived(), right.type().value());
    // }

    /**
     * @ingroup elementwise
     * @brief   Pointwise max.
     */
    template<class Left, class Right>
    UTOPIA_INLINE_FUNCTION DeviceBinary<Left, Right, Max> max(const DeviceExpression<Left> &left, const DeviceExpression<Right> &right) {
        return DeviceBinary<Left, Right, Max>(left.derived(), right.derived());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise max.
     */
    // template<class Left, class Right, int Order>
    // UTOPIA_INLINE_FUNCTION DeviceBinary<Left, DeviceNumber<Right>, Max> max(const DeviceExpression<Left> &left, const Factory<Values<Right>, Order> &right) {
    //     return DeviceBinary<Left, DeviceNumber<Right>, Max>(left.derived(), right.type().value());
    // }

    /**
     * @ingroup elementwise
     * @brief   Pointwise max.
     */
    // template<class Left, class Right, int Order>
    // UTOPIA_INLINE_FUNCTION DeviceBinary<Left, DeviceNumber<Right>, Max> max(const Factory<Values<Right>, Order> &right, const DeviceExpression<Left> &left) {
    //     return DeviceBinary<Left, DeviceNumber<Right>, Max>(left.derived(), right.type().value());
    // }

    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceTranspose<Derived> transpose(const DeviceExpression<Derived> &expr)
    {
        return DeviceTranspose<Derived>(expr.derived());
    }

    template<class Derived>
    UTOPIA_INLINE_FUNCTION UTOPIA_STORE_CONST(Derived) transpose(const DeviceTranspose<Derived> &expr)
    {
        return expr.expr();
    }

    template<class Expr>
    UTOPIA_INLINE_FUNCTION typename Traits<Expr>::Scalar det(
        const DeviceExpression<Expr> &expr)
    {
        return DeviceDeterminant<Expr>::apply(expr.derived());
    }

    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceInverse<Derived> inv(const DeviceExpression<Derived> &expr)
    {
        return DeviceInverse<Derived>(expr.derived());
    }

    template<class Derived>
    UTOPIA_INLINE_FUNCTION UTOPIA_STORE_CONST(Derived) inv(const DeviceInverse<Derived> &expr)
    {
        return expr.expr();
    }


    /////////////////////////////////////

    template<class Expr, class Result>
    UTOPIA_INLINE_FUNCTION void eig(
        const DeviceExpression<Expr> &expr,
        Result &eigen_values
        )
    {
        DeviceEigenValues<Expr>::apply(expr.derived(), eigen_values);
    }

    template<class Expr, class Values, class Vectors>
    UTOPIA_INLINE_FUNCTION void eig(
        const DeviceExpression<Expr> &expr,
        Values &eigen_values,
        Vectors &eigen_vectors
        )
    {
        DeviceEigenValues<Expr>::apply(expr.derived(), eigen_values);
        DeviceEigenVectors<Expr>::apply(expr.derived(), eigen_values, eigen_vectors);
    }


    template<class Expr, class Result>
    UTOPIA_INLINE_FUNCTION void sv(
        const DeviceExpression<Expr> &expr,
        Result &singular_values
        )
    {
        DeviceSingularValues<Expr>::apply(expr.derived(), singular_values);
    }

    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceDiag<Derived> diag(const DeviceExpression<Derived> &expr)
    {
        return DeviceDiag<Derived>(expr.derived());
    }

    template<class Derived>
    UTOPIA_INLINE_FUNCTION typename Traits<Derived>::Scalar sum(const DeviceExpression<Derived> &expr)
    {
        return DeviceReduce<Derived, Plus, Traits<Derived>::Order>::apply(expr.derived(), 0.0);
    }

    template<class Derived>
    UTOPIA_INLINE_FUNCTION DeviceTensorReduce<Derived, Plus, 1> row_sum(const DeviceExpression<Derived> &expr)
    {
        return DeviceTensorReduce<Derived, Plus, 1>(expr.derived());
    }
}

#endif //UTOPIA_DEVICE_OPERATIONS_HPP
