#ifndef UTOPIA_VIEW_OPERATIONS_HPP
#define UTOPIA_VIEW_OPERATIONS_HPP


#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Operators.hpp"
#include "utopia_DeviceUnary.hpp"
#include "utopia_DeviceNegate.hpp"
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
    DeviceBinary<Multiply<Tensor<L, 2>, Tensor<R, 1>>, Tensor<Left, 1>, Plus> operator+(const Tensor<Left, 1> &left, const Multiply<Tensor<L, 2>, Tensor<R, 1>> &right) {
        return DeviceBinary<Multiply<Tensor<L, 2>, Tensor<R, 1>>, Tensor<Left, 1>, Plus>(right.derived(), left.derived());
    }

    template<class Left, class Right>
    Multiply<Left, Right> operator*(const DeviceExpression<Left> &left, const DeviceExpression<Right> &right) {
        return Multiply<Left, Right>(left.derived(), right.derived());
    }

    template<class Left>
    DeviceBinary<Left, DeviceNumber<typename Left::Scalar>, Multiplies> operator*(const DeviceExpression<Left> &left,
                                                                       const typename Left::Scalar &right) {
        return DeviceBinary<Left, DeviceNumber<typename Left::Scalar>, Multiplies>(left.derived(), right);
    }

    template<class Right>
    DeviceBinary<DeviceNumber<typename Right::Scalar>, Right, Multiplies> operator*(const typename Right::Scalar &left,
                                                                         const DeviceExpression<Right> &right) {
        return DeviceBinary<DeviceNumber<typename Right::Scalar>, Right, Multiplies>(left, right.derived());
    }

    template<class Left, class Right, class Operation>
    DeviceBinary<Left, DeviceReduce<Right, Operation>, Multiplies> operator*(const DeviceExpression<Left> &left,
                                                                 const DeviceReduce<Right, Operation> &right) {
        return DeviceBinary<Left, DeviceReduce<Right, Operation>, Multiplies>(left.derived(), right);
    }

    template<class Left, class Right, int Order>
    DeviceBinary<Left, DeviceNorm<Right, Order>, Multiplies> operator*(const DeviceExpression<Left> &left,
                                                           const DeviceNorm<Right, Order> &right) {
        return DeviceBinary<Left, DeviceNorm<Right, Order>, Multiplies>(left.derived(), right);
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

    //Counteract the reciprocal capture of the / operator
    template<class Expr, class Operation>
    DeviceBinary<DeviceNumber<typename Expr::Scalar>,
           DeviceReduce<Expr, Operation>, Divides> operator / (const typename Expr::Scalar &left,
                                                         const DeviceReduce<Expr, Operation> &right)
    {
        return DeviceBinary<DeviceNumber<typename Expr::Scalar>,
                      DeviceReduce<Expr, Operation>, Divides>(left, right);
    }

    template<class Derived, class Expr, class Operation>
    DeviceBinary<Derived,
           DeviceReduce<Expr, Operation>, Divides> operator / (const DeviceExpression<Derived> &left,
                                                         const DeviceReduce<Expr, Operation> &right)
    {
        return DeviceBinary<Derived,
                      DeviceReduce<Expr, Operation>, Divides>(left.derived(), right);
    }

    template<class Expr, int Type>
    DeviceBinary<DeviceNumber<typename Expr::Scalar>,
           DeviceNorm<Expr, Type>, Divides> operator/(const typename Expr::Scalar &left,
                                                const DeviceNorm<Expr, Type> &right)
    {
        return DeviceBinary<DeviceNumber<typename Expr::Scalar>,
                      DeviceNorm<Expr, Type>, Divides>(left, right);
    }

    template<class LExpr, class LOp, class RExpr, class ROp>
    DeviceBinary<DeviceReduce<LExpr, LOp>,
           DeviceReduce<RExpr, ROp>, Divides> operator/(const DeviceReduce<LExpr, LOp> &left,
                                                  const DeviceReduce<RExpr, ROp> &right)
    {
        return DeviceBinary<DeviceReduce<LExpr, LOp>,
                      DeviceReduce<RExpr, ROp>, Divides>(left, right);
    }

    template<class Left, class Right>
    using EWApproxEqual = utopia::DeviceBinary<Left, Right, ApproxEqual>;

    template<class Left, class Right>
    using DeviceReduceApproxEqual = utopia::DeviceReduce<EWApproxEqual<Left, Right>, And>;

    template<class Left, class Right>
    Boolean<DeviceReduceApproxEqual<Left, Right>> approxeq(const DeviceExpression<Left> &left,
                                                     const DeviceExpression<Right> &right,
                                                     const typename Right::Scalar tol = 1e-6) {
        
        typedef utopia::EWApproxEqual<Left, Right> BinOp;
        return DeviceReduceApproxEqual<Left, Right>(BinOp(left.derived(), right.derived(), ApproxEqual(tol)));
    }

    /**     @defgroup   queries Structural and numerical queries
     *       @ingroup    algebra
    */

     /**
      * @ingroup queries
      * @brief   Checks if 2 variables are same up to requested tolerance.
      *
      * @param[in]  left   The left.
      * @param[in]  right  The right.
      * @param[in]  tol    The tolerance.
      */
    inline bool approxeq(const double left, const double right, const double tol = 10. * std::numeric_limits<double>::epsilon())
    {
        return std::abs(left-right) < tol;
    }

   inline bool approxeq(const float left, const float right, const double tol = 10.f * std::numeric_limits<float>::epsilon())
    {
        return std::abs(left-right) < tol;
    }


    /**     @defgroup   elementwise Element-wise
     *       @ingroup    algebra
    */



    /**
     * @ingroup elementwise
     * @brief   Pointwise multiplication.
     */
    template<class Left, class Right>
    inline DeviceBinary<Left, Right, EMultiplies> e_mul(const DeviceExpression<Left> &left, const DeviceExpression<Right> &right) {
        return DeviceBinary<Left, Right, EMultiplies>(left.derived(),
                                                right.derived());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise min.
     */
    template<class Left, class Right>
    inline DeviceBinary<Left, Right, Min> min(const DeviceExpression<Left> &left, const DeviceExpression<Right> &right) {
        return DeviceBinary<Left, Right, Min>(left.derived(), right.derived());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise min.
     */
    template<class Left, class Right, int Order>
    inline DeviceBinary<Left, DeviceNumber<Right>, Min> min(const DeviceExpression<Left> &left, const Factory<Values<Right>, Order> &right) {
        return DeviceBinary<Left, DeviceNumber<Right>, Min>(left.derived(), right.type().value());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise min.
     */
    template<class Left, class Right, int Order>
    inline DeviceBinary<Left, DeviceNumber<Right>, Min> min(const Factory<Values<Right>, Order> &right, const DeviceExpression<Left> &left) {
        return DeviceBinary<Left, DeviceNumber<Right>, Min>(left.derived(), right.type().value());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise max.
     */
    template<class Left, class Right>
    inline DeviceBinary<Left, Right, Max> max(const DeviceExpression<Left> &left, const DeviceExpression<Right> &right) {
        return DeviceBinary<Left, Right, Max>(left.derived(), right.derived());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise max.
     */
    template<class Left, class Right, int Order>
    inline DeviceBinary<Left, DeviceNumber<Right>, Max> max(const DeviceExpression<Left> &left, const Factory<Values<Right>, Order> &right) {
        return DeviceBinary<Left, DeviceNumber<Right>, Max>(left.derived(), right.type().value());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise max.
     */
    template<class Left, class Right, int Order>
    inline DeviceBinary<Left, DeviceNumber<Right>, Max> max(const Factory<Values<Right>, Order> &right, const DeviceExpression<Left> &left) {
        return DeviceBinary<Left, DeviceNumber<Right>, Max>(left.derived(), right.type().value());
    }
}

#endif //UTOPIA_VIEW_OPERATIONS_HPP
