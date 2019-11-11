#ifndef UTOPIA_VIEW_OPERATIONS_HPP
#define UTOPIA_VIEW_OPERATIONS_HPP


#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Operators.hpp"
#include "utopia_ViewUnary.hpp"
#include "utopia_ViewNegate.hpp"
#include "utopia_ViewBinary.hpp"
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
    ViewBinary<Expr, DeviceNumber<double>, Minus> operator-(const ViewExpression<Expr> &left, const double &right)
    {
        return ViewBinary<Expr, DeviceNumber<double>, Minus>(left, right);
    }

    template<class Expr>
    double operator-(const DeviceNumber<double> &left, const DeviceNumber<double> &right)
    {
        return static_cast<double>(left) - static_cast<double>(right);
    }

    template<class Derived>
    ViewNegate<Derived> operator-(const ViewExpression<Derived> &expr) {
        return ViewNegate<Derived>(expr.derived());
    }

    template<class Left, class Right>
    ViewBinary<Left, Right, Minus> operator-(const ViewExpression<Left> &left, const ViewExpression<Right> &right) {
        return ViewBinary<Left, Right, Minus>(left.derived(), right.derived());
    }


    template<class Left, class Right>
    ViewBinary<Left, Right, Plus> operator+(const ViewExpression<Left> &left, const ViewExpression<Right> &right) {
        return ViewBinary<Left, Right, Plus>(left.derived(), right.derived());
    }

    //Switch left with right 
    template<class Left, class L, class R>
    ViewBinary<Multiply<Tensor<L, 2>, Tensor<R, 1>>, Tensor<Left, 1>, Plus> operator+(const Tensor<Left, 1> &left, const Multiply<Tensor<L, 2>, Tensor<R, 1>> &right) {
        return ViewBinary<Multiply<Tensor<L, 2>, Tensor<R, 1>>, Tensor<Left, 1>, Plus>(right.derived(), left.derived());
    }

    template<class Left, class Right>
    Multiply<Left, Right> operator*(const ViewExpression<Left> &left, const ViewExpression<Right> &right) {
        return Multiply<Left, Right>(left.derived(), right.derived());
    }

    template<class Left>
    ViewBinary<Left, DeviceNumber<typename Left::Scalar>, Multiplies> operator*(const ViewExpression<Left> &left,
                                                                       const typename Left::Scalar &right) {
        return ViewBinary<Left, DeviceNumber<typename Left::Scalar>, Multiplies>(left.derived(), right);
    }

    template<class Right>
    ViewBinary<DeviceNumber<typename Right::Scalar>, Right, Multiplies> operator*(const typename Right::Scalar &left,
                                                                         const ViewExpression<Right> &right) {
        return ViewBinary<DeviceNumber<typename Right::Scalar>, Right, Multiplies>(left, right.derived());
    }

    template<class Left, class Right, class Operation>
    ViewBinary<Left, ViewReduce<Right, Operation>, Multiplies> operator*(const ViewExpression<Left> &left,
                                                                 const ViewReduce<Right, Operation> &right) {
        return ViewBinary<Left, ViewReduce<Right, Operation>, Multiplies>(left.derived(), right);
    }

    template<class Left, class Right, int Order>
    ViewBinary<Left, ViewNorm<Right, Order>, Multiplies> operator*(const ViewExpression<Left> &left,
                                                           const ViewNorm<Right, Order> &right) {
        return ViewBinary<Left, ViewNorm<Right, Order>, Multiplies>(left.derived(), right);
    }


    template<class Derived, class Operation>
    ViewUnary<Derived, Operation> transform(const ViewExpression<Derived> &expr, const Operation operation = Operation()) {
        return ViewUnary<Derived, Operation>(expr, operation);
    }

    /**     @defgroup   transforms Transforms
     *       @ingroup    algebra
    */


    /**
     * @ingroup transforms
     * @brief   \f$ x_i^2  \f$.
     */
    template<class Derived>
    auto pow2(const ViewExpression<Derived> &expr) -> decltype(transform(expr, Pow2())) {
        return transform(expr, Pow2());
    }
    template<class Derived>
    auto power(const ViewExpression<Derived> &expr, const double &a) -> decltype(transform(expr, Pow(a))) {
        return transform(expr, Pow(a));
    }

    template<class Derived>
    auto logn(const ViewExpression<Derived> &expr) -> decltype(transform(expr, Log())) {
        return transform(expr, Log());
    }

    template<class Derived>
    auto exp(const ViewExpression<Derived> &expr) -> decltype(transform(expr, Exp())) {
        return transform(expr, Exp());
    }

    template<class Derived>
    auto cos(const ViewExpression<Derived> &expr) -> decltype(transform(expr, Cos())) {
        return transform(expr, Cos());
    }

    template<class Derived>
    auto sin(const ViewExpression<Derived> &expr) -> decltype(transform(expr, Sin())) {
        return transform(expr, Sin());
    }


    /**
     * @ingroup transforms
     * @brief   \f$ | x_i |  \f$.
     */
    template<class Derived>
    auto abs(const ViewExpression<Derived> &expr) -> decltype(transform(expr, Abs())) {
        return transform(expr, Abs());
    }

    /**
     * @ingroup transforms
     * @brief   \f$ \sqrt(x_i)  \f$.
     */
    template<class Derived>
    auto sqrt(const ViewExpression<Derived> &expr) -> decltype(transform(expr, Sqrt())) {
        return transform(expr, Sqrt());
    }

    /**
     * @ingroup transforms
     * @brief   \f$ \frac 1 {x_i}  \f$.
     */
    template<class Derived>
    auto operator / (const typename Derived::Scalar &left,
                     const ViewExpression<Derived> &expr) -> decltype( transform(expr, Reciprocal<typename Derived::Scalar>(left)) ) {
        return transform(expr, Reciprocal<typename Derived::Scalar>(left));
    }

    template<class LDerived, class RDerived>
    ViewBinary<LDerived, RDerived, Divides> operator / (const ViewExpression<LDerived> &left,
                                                    const ViewExpression<RDerived> &right)
    {
        return ViewBinary<LDerived, RDerived, Divides>(left.derived(), right.derived());
    }

    //Counteract the reciprocal capture of the / operator
    template<class Expr, class Operation>
    ViewBinary<DeviceNumber<typename Expr::Scalar>,
           ViewReduce<Expr, Operation>, Divides> operator / (const typename Expr::Scalar &left,
                                                         const ViewReduce<Expr, Operation> &right)
    {
        return ViewBinary<DeviceNumber<typename Expr::Scalar>,
                      ViewReduce<Expr, Operation>, Divides>(left, right);
    }

    template<class Derived, class Expr, class Operation>
    ViewBinary<Derived,
           ViewReduce<Expr, Operation>, Divides> operator / (const ViewExpression<Derived> &left,
                                                         const ViewReduce<Expr, Operation> &right)
    {
        return ViewBinary<Derived,
                      ViewReduce<Expr, Operation>, Divides>(left.derived(), right);
    }

    template<class Expr, int Type>
    ViewBinary<DeviceNumber<typename Expr::Scalar>,
           ViewNorm<Expr, Type>, Divides> operator/(const typename Expr::Scalar &left,
                                                const ViewNorm<Expr, Type> &right)
    {
        return ViewBinary<DeviceNumber<typename Expr::Scalar>,
                      ViewNorm<Expr, Type>, Divides>(left, right);
    }

    template<class LExpr, class LOp, class RExpr, class ROp>
    ViewBinary<ViewReduce<LExpr, LOp>,
           ViewReduce<RExpr, ROp>, Divides> operator/(const ViewReduce<LExpr, LOp> &left,
                                                  const ViewReduce<RExpr, ROp> &right)
    {
        return ViewBinary<ViewReduce<LExpr, LOp>,
                      ViewReduce<RExpr, ROp>, Divides>(left, right);
    }

    template<class Left, class Right>
    using EWApproxEqual = utopia::ViewBinary<Left, Right, ApproxEqual>;

    template<class Left, class Right>
    using ViewReduceApproxEqual = utopia::ViewReduce<EWApproxEqual<Left, Right>, And>;

    template<class Left, class Right>
    Boolean<ViewReduceApproxEqual<Left, Right>> approxeq(const ViewExpression<Left> &left,
                                                     const ViewExpression<Right> &right,
                                                     const typename Right::Scalar tol = 1e-6) {
        
        typedef utopia::EWApproxEqual<Left, Right> BinOp;
        return ViewReduceApproxEqual<Left, Right>(BinOp(left.derived(), right.derived(), ApproxEqual(tol)));
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
    inline ViewBinary<Left, Right, EMultiplies> e_mul(const ViewExpression<Left> &left, const ViewExpression<Right> &right) {
        return ViewBinary<Left, Right, EMultiplies>(left.derived(),
                                                right.derived());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise min.
     */
    template<class Left, class Right>
    inline ViewBinary<Left, Right, Min> min(const ViewExpression<Left> &left, const ViewExpression<Right> &right) {
        return ViewBinary<Left, Right, Min>(left.derived(), right.derived());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise min.
     */
    template<class Left, class Right, int Order>
    inline ViewBinary<Left, DeviceNumber<Right>, Min> min(const ViewExpression<Left> &left, const Factory<Values<Right>, Order> &right) {
        return ViewBinary<Left, DeviceNumber<Right>, Min>(left.derived(), right.type().value());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise min.
     */
    template<class Left, class Right, int Order>
    inline ViewBinary<Left, DeviceNumber<Right>, Min> min(const Factory<Values<Right>, Order> &right, const ViewExpression<Left> &left) {
        return ViewBinary<Left, DeviceNumber<Right>, Min>(left.derived(), right.type().value());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise max.
     */
    template<class Left, class Right>
    inline ViewBinary<Left, Right, Max> max(const ViewExpression<Left> &left, const ViewExpression<Right> &right) {
        return ViewBinary<Left, Right, Max>(left.derived(), right.derived());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise max.
     */
    template<class Left, class Right, int Order>
    inline ViewBinary<Left, DeviceNumber<Right>, Max> max(const ViewExpression<Left> &left, const Factory<Values<Right>, Order> &right) {
        return ViewBinary<Left, DeviceNumber<Right>, Max>(left.derived(), right.type().value());
    }

    /**
     * @ingroup elementwise
     * @brief   Pointwise max.
     */
    template<class Left, class Right, int Order>
    inline ViewBinary<Left, DeviceNumber<Right>, Max> max(const Factory<Values<Right>, Order> &right, const ViewExpression<Left> &left) {
        return ViewBinary<Left, DeviceNumber<Right>, Max>(left.derived(), right.type().value());
    }
}

#endif //UTOPIA_VIEW_OPERATIONS_HPP
