#ifndef UTOPIA_HAS_MAT_MAT_MUL_HPP
#define UTOPIA_HAS_MAT_MAT_MUL_HPP

#include "utopia_Expressions.hpp"

namespace utopia {
    // template<class Expr>
    // class HasMatMatMul {
    // public:
    // 	enum { value = 0 };
    // };

    // template<class InnerExpr, class Operation>
    // class HasMatMatMul< Unary<InnerExpr, Operation> > {
    // public:
    // 	enum { value = HasMatMatMul<InnerExpr>::value };
    // };

    // template<class InnerExpr, class Operation>
    // class HasMatMatMul< Reduce<InnerExpr, Operation> > {
    // public:
    // 	enum { value = HasMatMatMul<InnerExpr>::value };
    // };

    // template<class Left, class Right, class Operation>
    // class HasMatMatMul< Binary<Left, Right, Operation> > {
    // public:
    // 	enum { value = HasMatMatMul<Left>::value || HasMatMatMul<Right>::value };
    // };

    // template<class Left, class Right>
    // class HasMatMatMul< Multiply<Left, Right> > {
    // public:
    // 	enum { value = HasMatMatMul<Left>::value || HasMatMatMul<Right>::value || ( Left::Order == 2 && Right::Order >= 1) };
    // };

    // template<class InnerExpr>
    // class HasMatMatMul< Transposed<InnerExpr> > {
    // public:
    // 	enum { value = HasMatMatMul<InnerExpr>::value };
    // };

    // template<class Derived>
    // inline constexpr bool has_mat_mat_mul(const Expression<Derived> &)
    // {
    // 	return HasMatMatMul<Derived>::value;
    // }
}

#endif //UTOPIA_HAS_MAT_MAT_MUL_HPP

