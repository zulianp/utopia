#ifndef UTOPIA_FE_IS_SUBTREE_HPP
#define UTOPIA_FE_IS_SUBTREE_HPP

#include "utopia_Any.hpp"
#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_IsSubTree.hpp"
#include "utopia_Operators.hpp"

namespace utopia {

    // trial function
    template <class Space>
    class IsSubTree<TrialFunction<Space>, TrialFunction<Space>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Space>
    class IsSubTree<Expr, TrialFunction<Space>> {
    public:
        static const int value = 0;
    };

    template <class Space>
    class IsSubTree<TrialFunction<utopia::Any>, TrialFunction<Space>> {
    public:
        static const int value = 1;
    };

    // test function
    template <class Space>
    class IsSubTree<TestFunction<Space>, TestFunction<Space>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Space>
    class IsSubTree<Expr, TestFunction<Space>> {
    public:
        static const int value = 0;
    };

    template <class Space>
    class IsSubTree<TestFunction<utopia::Any>, TestFunction<Space>> {
    public:
        static const int value = 1;
    };

    // gradient
    template <class Inner>
    class IsSubTree<Gradient<Inner>, Gradient<Inner>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner>
    class IsSubTree<Expr, Gradient<Inner>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // test gradient
    template <class Space>
    class IsSubTree<Gradient<TestFunction<utopia::Any>>, Gradient<TestFunction<Space>>> {
    public:
        static const int value = 1;
    };

    // trial gradient
    template <class Space>
    class IsSubTree<Gradient<TrialFunction<utopia::Any>>, Gradient<TrialFunction<Space>>> {
    public:
        static const int value = 1;
    };

    // divergence
    template <class Inner>
    class IsSubTree<Divergence<Inner>, Divergence<Inner>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner>
    class IsSubTree<Expr, Divergence<Inner>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // test divergence
    template <class Space>
    class IsSubTree<Divergence<TestFunction<utopia::Any>>, Divergence<TestFunction<Space>>> {
    public:
        static const int value = 1;
    };

    // trial divergence
    template <class Space>
    class IsSubTree<Divergence<TrialFunction<utopia::Any>>, Divergence<TrialFunction<Space>>> {
    public:
        static const int value = 1;
    };

    // curl
    template <class Inner>
    class IsSubTree<Curl<Inner>, Curl<Inner>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner>
    class IsSubTree<Expr, Curl<Inner>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // test curl
    template <class Space>
    class IsSubTree<Curl<TestFunction<utopia::Any>>, Curl<TestFunction<Space>>> {
    public:
        static const int value = 1;
    };

    // trial curl
    template <class Space>
    class IsSubTree<Curl<TrialFunction<utopia::Any>>, Curl<TrialFunction<Space>>> {
    public:
        static const int value = 1;
    };

    // time-derivative
    template <class Inner, int Order>
    class IsSubTree<TimeDerivative<Inner, Order>, TimeDerivative<Inner, Order>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner, int Order>
    class IsSubTree<Expr, TimeDerivative<Inner, Order>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    template <class Inner, int Order>
    class IsSubTree<TimeDerivative<utopia::Any, Order>, TimeDerivative<Inner, Order>> {
    public:
        static const int value = 1;
    };

    // integral
    template <class Inner>
    class IsSubTree<Integral<Inner>, Integral<Inner>> {
    public:
        static const int value = 1;
    };

    template <class Inner>
    class IsSubTree<Integral<utopia::Any>, Integral<Inner>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class Inner>
    class IsSubTree<Expr, Integral<Inner>> {
    public:
        static const int value = IsSubTree<Expr, Inner>::value;
    };

    // equality
    template <class Expr, class Left, class Right>
    class IsSubTree<Expr, Equality<Left, Right>> {
    public:
        static const int value = IsSubTree<Expr, Left>::value || IsSubTree<Expr, Right>::value;
    };

    template <class Left, class Right>
    class IsSubTree<Equality<Left, Right>, Equality<Left, Right>> {
    public:
        static const int value = 1;
    };

    template <class Left, class Right>
    class IsSubTree<Equality<utopia::Any, utopia::Any>, Equality<Left, Right>> {
    public:
        static const int value = 1;
    };

    template <class Expr, class First, class... Rest>
    class IsSubTree<Expr, Equations<First, Rest...>> {
    public:
        static const int first_value = IsSubTree<Expr, First>::value;
        static const int rest_value = IsSubTree<Expr, Equations<Rest...>>::value;
        static const int value = (first_value > rest_value) ? first_value : rest_value;
    };

    template <class Expr, class Eq>
    class IsSubTree<Expr, Equations<Eq>> {
    public:
        static const int value = IsSubTree<Expr, Eq>::value;
    };

    // interpolate
    template <class T1, class T2>
    class IsSubTree<Interpolate<utopia::Any, utopia::Any>, Interpolate<T1, T2>> {
    public:
        static const int value = 1;
    };

    // to be moved in utopia

    template <class Left, class Right>
    class IsSubTree<Reduce<Binary<Any, Any, EMultiplies>, Plus>, Reduce<Binary<Left, Right, EMultiplies>, Plus>> {
    public:
        static const int value = 1;
    };

    template <class Left, class Right>
    class IsSubTree<Binary<Any, Any, Multiplies>, Binary<Left, Right, Multiplies>> {
    public:
        static const int value = 1;
    };

    /// rest of specializations

    // filter
    template <class Inner, class Out, int Order, class Fun>
    class IsSubTree<Filter<Inner, Out, Order, Fun>, Filter<Inner, Out, Order, Fun>> {
    public:
        static const int value = 1;
    };

    template <class Inner, class Out, int Order, class Fun>
    class IsSubTree<Filter<utopia::Any, Out, Order, Fun>, Filter<Inner, Out, Order, Fun>> {
    public:
        static const int value = 1;
    };

    template <class Inner, class Out, int Order, class Fun>
    class IsSubTree<Inner, Filter<Inner, Out, Order, Fun>> {
    public:
        static const int value = IsSubTree<Inner, Inner>::value;
    };
}  // namespace utopia

#endif  // UTOPIA_FE_IS_SUBTREE_HPP
