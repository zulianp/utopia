//
// Created by Patrick Zulian on 15/05/15.
//

#ifndef utopia_utopia_OPERATIONS_HPP
#define utopia_utopia_OPERATIONS_HPP


#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Operators.hpp"
#include "utopia_Unary.hpp"
#include "utopia_Negate.hpp"
#include "utopia_Binary.hpp"
#include "utopia_Literal.hpp"
#include "utopia_Boolean.hpp"

#include <limits>
#include <cmath>

namespace utopia {
    template<class Expr>
    Binary<Expr, Number<double>, Minus> operator-(const Expression<Expr> &left, const double &right)
    {
        return Binary<Expr, Number<double>, Minus>(left, right); 
    }

    template<class Expr>
    double operator-(const Number<double> &left, const Number<double> &right)
    {
        return static_cast<double>(left) - static_cast<double>(right);
    }

    template<class Derived>
    Negate<Derived> operator-(const Expression<Derived> &expr) {
        return Negate<Derived>(expr.derived());
    }

    template<class Left, class Right>
    Binary<Left, Right, Minus> operator-(const Expression<Left> &left, const Expression<Right> &right) {
        return Binary<Left, Right, Minus>(left.derived(), right.derived());
    }


    template<class Left, class Right>
    Binary<Left, Right, Plus> operator+(const Expression<Left> &left, const Expression<Right> &right) {
        return Binary<Left, Right, Plus>(left.derived(), right.derived());
    }

    template<class Left, class Right>
    Multiply<Left, Right> operator*(const Expression<Left> &left, const Expression<Right> &right) {
        return Multiply<Left, Right>(left.derived(), right.derived());
    }

    template<class Left>
    Binary<Left, Number<typename Left::Scalar>, Multiplies> operator*(const Expression<Left> &left,
                                                                       const typename Left::Scalar &right) {
        return Binary<Left, Number<typename Left::Scalar>, Multiplies>(left.derived(), right);
    }

    template<class Right>
    Binary<Number<typename Right::Scalar>, Right, Multiplies> operator*(const typename Right::Scalar &left,
                                                                         const Expression<Right> &right) {
        return Binary<Number<typename Right::Scalar>, Right, Multiplies>(left, right.derived());
    }

    template<class Left, class Right, class Operation>
    Binary<Left, Reduce<Right, Operation>, Multiplies> operator*(const Expression<Left> &left,
                                                                 const Reduce<Right, Operation> &right) {
        return Binary<Left, Reduce<Right, Operation>, Multiplies>(left.derived(), right);
    }

    template<class Left, class Right, int Order>
    Binary<Left, Norm<Right, Order>, Multiplies> operator*(const Expression<Left> &left,
                                                           const Norm<Right, Order> &right) {
        return Binary<Left, Norm<Right, Order>, Multiplies>(left.derived(), right);
    }


    template<class Derived, class Operation>
    Unary<Derived, Operation> transform(const Expression<Derived> &expr, const Operation operation = Operation()) {
        return Unary<Derived, Operation>(expr, operation);
    }

    /**     @defgroup   transforms Transforms
     *       @ingroup    algebra
    */


    /**
     * @ingroup transforms
     * @brief   \f$ x_i^2  \f$.
     */
    template<class Derived>
    auto pow2(const Expression<Derived> &expr) -> decltype(transform(expr, Pow2())) {
        return transform(expr, Pow2());
    }
    template<class Derived>
    auto power(const Expression<Derived> &expr, const double &a) -> decltype(transform(expr, Pow(a))) {
        return transform(expr, Pow(a));
    }

    template<class Derived>
    auto logn(const Expression<Derived> &expr) -> decltype(transform(expr, Log())) {
        return transform(expr, Log());
    }

    template<class Derived>
    auto exp(const Expression<Derived> &expr) -> decltype(transform(expr, Exp())) {
        return transform(expr, Exp());
    }

    template<class Derived>
    auto cos(const Expression<Derived> &expr) -> decltype(transform(expr, Cos())) {
        return transform(expr, Cos());
    }

    template<class Derived>
    auto sin(const Expression<Derived> &expr) -> decltype(transform(expr, Sin())) {
        return transform(expr, Sin());
    }


    /**
     * @ingroup transforms
     * @brief   \f$ | x_i |  \f$.
     */
    template<class Derived>
    auto abs(const Expression<Derived> &expr) -> decltype(transform(expr, Abs())) {
        return transform(expr, Abs());
    }

    /**
     * @ingroup transforms
     * @brief   \f$ \sqrt(x_i)  \f$.
     */
    template<class Derived>
    auto sqrt(const Expression<Derived> &expr) -> decltype(transform(expr, Sqrt())) {
        return transform(expr, Sqrt());
    }

    /**
     * @ingroup transforms
     * @brief   \f$ \frac 1 {x_i}  \f$.
     */
    template<class Derived>
    auto operator / (const typename Derived::Scalar &left, 
                     const Expression<Derived> &expr) -> decltype( transform(expr, Reciprocal<typename Derived::Scalar>(left)) ) {
        return transform(expr, Reciprocal<typename Derived::Scalar>(left));
    }

    template<class LDerived, class RDerived>
    Binary<LDerived, RDerived, Divides> operator / (const Expression<LDerived> &left,
                                                    const Expression<RDerived> &right) 
    {
        return Binary<LDerived, RDerived, Divides>(left.derived(), right.derived());
    }

    //Counteract the reciprocal capture of the / operator
    template<class Expr, class Operation>
    Binary<Number<typename Expr::Scalar>, 
           Reduce<Expr, Operation>, Divides> operator / (const typename Expr::Scalar &left, 
                                                         const Reduce<Expr, Operation> &right) 
    {
        return Binary<Number<typename Expr::Scalar>, 
                      Reduce<Expr, Operation>, Divides>(left, right);
    }

    template<class Derived, class Expr, class Operation>
    Binary<Derived, 
           Reduce<Expr, Operation>, Divides> operator / (const Expression<Derived> &left, 
                                                         const Reduce<Expr, Operation> &right) 
    {
        return Binary<Derived, 
                      Reduce<Expr, Operation>, Divides>(left.derived(), right);
    }

    template<class Expr, int Type>
    Binary<Number<typename Expr::Scalar>, 
           Norm<Expr, Type>, Divides> operator/(const typename Expr::Scalar &left, 
                                                const Norm<Expr, Type> &right) 
    {
        return Binary<Number<typename Expr::Scalar>, 
                      Norm<Expr, Type>, Divides>(left, right);
    }

    template<class LExpr, class LOp, class RExpr, class ROp>
    Binary<Reduce<LExpr, LOp>, 
           Reduce<RExpr, ROp>, Divides> operator/(const Reduce<LExpr, LOp> &left,
                                                  const Reduce<RExpr, ROp> &right) 
    {
        return Binary<Reduce<LExpr, LOp>, 
                      Reduce<RExpr, ROp>, Divides>(left, right);
    }

    template<class Left, class Right>
    Boolean<Reduce<Binary<Left, Right, ApproxEqual>, And> > approxeq(const Expression<Left> &left,
                                                                     const Expression<Right> &right,
                                                           const typename Right::Scalar tol = 1e-6) {
        typedef utopia::Binary<Left, Right, ApproxEqual> BinOp;
        return Reduce<BinOp, And>(BinOp(left.derived(), right.derived(), ApproxEqual(tol)));
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


    
    /**
     * @ingroup tensor_products
     * @brief   Pointwise multiplication. 
     */
    template<class Left, class Right>
    inline Binary<Left, Right, EMultiplies> e_mul(const Expression<Left> &left, const Expression<Right> &right) {
        return Binary<Left, Right, EMultiplies>(left.derived(), 
                                                right.derived());
    }

    /**
     * @ingroup tensor_products
     * @brief   Pointwise min. 
     */
    template<class Left, class Right>
    inline Binary<Left, Right, Min> min(const Expression<Left> &left, const Expression<Right> &right) {
        return Binary<Left, Right, Min>(left.derived(), right.derived());
    }


    /**
     * @ingroup tensor_products
     * @brief   Pointwise max. 
     */
    template<class Left, class Right>
    inline Binary<Left, Right, Max> max(const Expression<Left> &left, const Expression<Right> &right) {
        return Binary<Left, Right, Max>(left.derived(), right.derived());
    }


}

#endif //utopia_utopia_OPERATIONS_HPP
