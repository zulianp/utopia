//
// Created by Patrick Zulian on 03/06/15.
//

#ifndef UTOPIA_UTOPIA_STOREAS_HPP
#define UTOPIA_UTOPIA_STOREAS_HPP


#define UTOPIA_STORE(_Expr) typename utopia::StoreAs<_Expr, _Expr::StoreAs, !UTOPIA_CONST>::Expr
#define UTOPIA_STORE_CONST(_Expr) typename utopia::StoreAs<_Expr, _Expr::StoreAs, UTOPIA_CONST>::Expr


namespace utopia {
    static const int UTOPIA_BY_VALUE = 1;
    static const int UTOPIA_BY_REFERENCE = 0;
    static const int UTOPIA_CONST = 1;


    template<typename _Expr, int ExprStorage = UTOPIA_BY_REFERENCE, bool Const = UTOPIA_CONST>
    struct StoreAs {
    };

    template<typename _Expr>
    struct StoreAs<_Expr, UTOPIA_BY_REFERENCE, !UTOPIA_CONST> {
        typedef _Expr &Expr;
    };

    template<typename _Expr>
    struct StoreAs<_Expr, UTOPIA_BY_VALUE, !UTOPIA_CONST> {
        typedef _Expr Expr;
    };

    template<typename _Expr>
    struct StoreAs<_Expr, UTOPIA_BY_REFERENCE, UTOPIA_CONST> {
        typedef const _Expr &Expr;
    };

    template<typename _Expr>
    struct StoreAs<_Expr, UTOPIA_BY_VALUE, UTOPIA_CONST> {
        typedef const _Expr Expr;
    };

    // static const int UTOPIA_DEFAULT_EXPRESSION_STORAGE = UTOPIA_BY_REFERENCE;
    static const int UTOPIA_DEFAULT_EXPRESSION_STORAGE = UTOPIA_BY_VALUE; ///much safer version of storing expressions (when subclassing Matrix add UTOPIA_BY_REFERENCE)

} /* express */

#endif //UTOPIA_UTOPIA_STOREAS_HPP
