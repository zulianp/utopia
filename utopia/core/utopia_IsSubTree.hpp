#ifndef UTOPIA_TREE_HAS_EXPRESSION_HPP
#define UTOPIA_TREE_HAS_EXPRESSION_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    // generic
    template <class Expr, class Tree>
    class IsSubTree {
    public:
        static const int value = 0;
    };

    template <class Expr>
    class IsSubTree<Expr, Expr> {
    public:
        static const int value = 1;
    };

    // multiply
    template <class Expr, class Left, class Right>
    class IsSubTree<Expr, Multiply<Left, Right>> {
    public:
        static const int value = IsSubTree<Expr, Left>::value || IsSubTree<Expr, Right>::value;
    };

    template <class Left, class Right>
    class IsSubTree<Multiply<Left, Right>, Multiply<Left, Right>> {
    public:
        static const int value = 1;
    };

    // binary
    template <class Left, class Right, class Op>
    class IsSubTree<Binary<Left, Right, Op>, Binary<Left, Right, Op>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Left, class Right, class Op>
    class IsSubTree<Expr, Binary<Left, Right, Op>> {
    public:
        static const int value = IsSubTree<Expr, Left>::value || IsSubTree<Expr, Right>::value;
    };

    // unary
    template <class Inner, class Op>
    class IsSubTree<Unary<Inner, Op>, Unary<Inner, Op>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner, class Op>
    class IsSubTree<Expr, Unary<Inner, Op>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // negate
    template <class Inner>
    class IsSubTree<Negate<Inner>, Negate<Inner>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner>
    class IsSubTree<Expr, Negate<Inner>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // reduce
    template <class Inner, class Op>
    class IsSubTree<Reduce<Inner, Op>, Reduce<Inner, Op>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner, class Op>
    class IsSubTree<Expr, Reduce<Inner, Op>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // tensor-reduce
    template <class Inner, class Op>
    class IsSubTree<TensorReduce<Inner, Op>, TensorReduce<Inner, Op>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner, class Op>
    class IsSubTree<Expr, TensorReduce<Inner, Op>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // transposed
    template <class Inner>
    class IsSubTree<Transposed<Inner>, Transposed<Inner>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner>
    class IsSubTree<Expr, Transposed<Inner>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // trace
    template <class Inner>
    class IsSubTree<Trace<Inner>, Trace<Inner>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner>
    class IsSubTree<Expr, Trace<Inner>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // inverse
    template <class Inner>
    class IsSubTree<Inverse<Inner>, Inverse<Inner>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner>
    class IsSubTree<Expr, Inverse<Inner>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // determinant
    template <class Inner>
    class IsSubTree<Determinant<Inner>, Determinant<Inner>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner>
    class IsSubTree<Expr, Determinant<Inner>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // diagonal
    template <class Inner>
    class IsSubTree<Diag<Inner>, Diag<Inner>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner>
    class IsSubTree<Expr, Diag<Inner>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // norm
    template <class Inner, int Type>
    class IsSubTree<Norm<Inner, Type>, Norm<Inner, Type>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner, int Type>
    class IsSubTree<Expr, Norm<Inner, Type>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // factory
    template <class Type, int Order>
    class IsSubTree<Factory<Type, Order>, Factory<Type, Order>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Type, int Order>
    class IsSubTree<Expr, Factory<Type, Order>> {
    public:
        static const int value = 0;
    };

    // wrapper
    template <class T, int Order>
    class IsSubTree<Tensor<T, Order>, Tensor<T, Order>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class T, int Order>
    class IsSubTree<Expr, Tensor<T, Order>> {
    public:
        static const int value = IsSubTree<Expr, T>::value;
    };

    // queries
    template <template <class...> class Expr>
    class ExprQuery {
    public:
        static const int value = 0;

        template <class... T2>
        inline constexpr static bool is_specialization(const Expr<T2...> &) {
            return true;
        }

        template <class Any>
        inline constexpr static bool is_specialization(const Any &) {
            return false;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_TREE_HAS_EXPRESSION_HPP
