#ifndef UTOPIA_DEVICE_OPERATIONS_HPP
#define UTOPIA_DEVICE_OPERATIONS_HPP


#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Operators.hpp"
#include "utopia_DeviceUnary.hpp"
#include "utopia_DeviceBinary.hpp"
#include "utopia_Literal.hpp"
#include "utopia_Boolean.hpp"
#include "utopia_Factory.hpp"

#ifdef WITH_TRILINOS
#include <Kokkos_ArithTraits.hpp>
#endif

#include <limits>
#include <algorithm>
#include <cmath>

namespace utopia {

    template<class Expr>
    DeviceBinary<Expr, DeviceNumber<double>, Minus> operator-(const DeviceExpression<Expr> &left, const double &right)
    {
        return DeviceBinary<Expr, DeviceNumber<double>, Minus>(left, right);
    }

    template<class Expr>
    double operator-(const DeviceNumber<double> &left, const DeviceNumber<double> &right)
    {
        return static_cast<double>(left) - static_cast<double>(right);
    }

    template<class Derived>
    DeviceNegate<Derived> operator-(const DeviceExpression<Derived> &expr) {
        return DeviceNegate<Derived>(expr.derived());
    }

    template<class Left, class Right>
    DeviceBinary<Left, Right, Minus> operator-(const DeviceExpression<Left> &left, const DeviceExpression<Right> &right) {
        return DeviceBinary<Left, Right, Minus>(left.derived(), right.derived());
    }

    template<class Left, class Right>
    DeviceBinary<Left, Right, Plus> operator+(const DeviceExpression<Left> &left, const DeviceExpression<Right> &right) {
        return DeviceBinary<Left, Right, Plus>(left.derived(), right.derived());
    }

    //Switch left with right 
    template<class Left, class L, class R>
    DeviceBinary<DeviceMultiply<TensorView<L, 2>, TensorView<R, 1>>, TensorView<Left, 1>, Plus> operator+(const TensorView<Left, 1> &left, const DeviceMultiply<TensorView<L, 2>, TensorView<R, 1>> &right) {
        return DeviceBinary<DeviceMultiply<TensorView<L, 2>, TensorView<R, 1>>, TensorView<Left, 1>, Plus>(right.derived(), left.derived());
    }

    template<class Left, class Right>
    DeviceMultiply<Left, Right> operator*(
        const DeviceExpression<Left> &left,
        const DeviceExpression<Right> &right
    ) {
        return DeviceMultiply<Left, Right>(left.derived(), right.derived());
    }

    template<class Left>
    DeviceBinary<
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
    DeviceBinary<
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
    DeviceUnary<Derived, Operation> transform(const DeviceExpression<Derived> &expr, const Operation operation = Operation()) {
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
    auto pow2(const DeviceExpression<Derived> &expr) -> decltype(transform(expr, Pow2())) {
        return transform(expr, Pow2());
    }
    template<class Derived>
    auto power(const DeviceExpression<Derived> &expr, const double &a) -> decltype(transform(expr, Pow(a))) {
        return transform(expr, Pow(a));
    }

    template<class Derived>
    auto logn(const DeviceExpression<Derived> &expr) -> decltype(transform(expr, Log())) {
        return transform(expr, Log());
    }

    template<class Derived>
    auto exp(const DeviceExpression<Derived> &expr) -> decltype(transform(expr, Exp())) {
        return transform(expr, Exp());
    }

    template<class Derived>
    auto cos(const DeviceExpression<Derived> &expr) -> decltype(transform(expr, Cos())) {
        return transform(expr, Cos());
    }

    template<class Derived>
    auto sin(const DeviceExpression<Derived> &expr) -> decltype(transform(expr, Sin())) {
        return transform(expr, Sin());
    }

    /**
     * @ingroup transforms
     * @brief   \f$ | x_i |  \f$.
     */
    template<class Derived>
    auto abs(const DeviceExpression<Derived> &expr) -> decltype(transform(expr, Abs())) {
        return transform(expr, Abs());
    }

    /**
     * @ingroup transforms
     * @brief   \f$ \sqrt(x_i)  \f$.
     */
    template<class Derived>
    auto sqrt(const DeviceExpression<Derived> &expr) -> decltype(transform(expr, Sqrt())) {
        return transform(expr, Sqrt());
    }

    /**
     * @ingroup transforms
     * @brief   \f$ \frac 1 {x_i}  \f$.
     */
    template<class Derived>
    auto operator / (const typename Derived::Scalar &left,
                     const DeviceExpression<Derived> &expr) -> decltype( transform(expr, Reciprocal<typename Derived::Scalar>(left)) ) {
        return transform(expr, Reciprocal<typename Derived::Scalar>(left));
    }

    template<class LDerived, class RDerived>
    DeviceBinary<LDerived, RDerived, Divides> operator / (const DeviceExpression<LDerived> &left,
                                                    const DeviceExpression<RDerived> &right)
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
    template<class Left, class Right, int Order>
    UTOPIA_INLINE_FUNCTION DeviceBinary<Left, DeviceNumber<Right>, Min> min(const DeviceExpression<Left> &left, const Factory<Values<Right>, Order> &right) {
        return DeviceBinary<Left, DeviceNumber<Right>, Min>(left.derived(), right.type().value());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise min.
     */
    template<class Left, class Right, int Order>
    UTOPIA_INLINE_FUNCTION DeviceBinary<Left, DeviceNumber<Right>, Min> min(const Factory<Values<Right>, Order> &right, const DeviceExpression<Left> &left) {
        return DeviceBinary<Left, DeviceNumber<Right>, Min>(left.derived(), right.type().value());
    }

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
    template<class Left, class Right, int Order>
    UTOPIA_INLINE_FUNCTION DeviceBinary<Left, DeviceNumber<Right>, Max> max(const DeviceExpression<Left> &left, const Factory<Values<Right>, Order> &right) {
        return DeviceBinary<Left, DeviceNumber<Right>, Max>(left.derived(), right.type().value());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise max.
     */
    template<class Left, class Right, int Order>
    UTOPIA_INLINE_FUNCTION DeviceBinary<Left, DeviceNumber<Right>, Max> max(const Factory<Values<Right>, Order> &right, const DeviceExpression<Left> &left) {
        return DeviceBinary<Left, DeviceNumber<Right>, Max>(left.derived(), right.type().value());
    }
}

#endif //UTOPIA_DEVICE_OPERATIONS_HPP
