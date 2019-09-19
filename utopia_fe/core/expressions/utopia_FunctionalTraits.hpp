#ifndef UTOPIA_FE_FUNCTIONAL_TRAITS_HPP
#define UTOPIA_FE_FUNCTIONAL_TRAITS_HPP

#include <iostream>
#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_Any.hpp"


namespace utopia {
    static const int CONSTANT_FUNCTION = -1;
    static const int POLYNOMIAL_FUNCTION = 0;
    static const int EXPONENTIAL_FUNCTION = 1;
    static const int UNDEFINED_FUNCTION = 2;

    template<class Expr, class AssemblyContext>
    class FunctionalTraits {
    public:
#ifdef NDEBUG
        inline constexpr static int type(const Expr &,  const AssemblyContext &) { return POLYNOMIAL_FUNCTION; }
        inline constexpr static int order(const Expr &, const AssemblyContext &) { return 0; }
#else
        inline static int type(const Expr &expr,  const AssemblyContext &)
        {
            std::cout << "Unhandled FunctionalTraits type for " <<  expr.get_class() << std::endl;
            return POLYNOMIAL_FUNCTION;
        }

        inline static int order(const Expr &expr, const AssemblyContext &)
        {
            std::cout << "Unhandled FunctionalTraits order for " <<  expr.get_class() << std::endl;
            return 0;
        }
#endif

    };

    template<class Left, class Right>
    class SelectFunctionSpace {
    public:
        typedef Left Type;
    };

    template<class Right>
    class SelectFunctionSpace<None, Right> {
    public:
        typedef Right Type;
    };

    template<class Left>
    class SelectFunctionSpace<Left, None> {
    public:
        typedef Left Type;
    };

    template<>
    class SelectFunctionSpace<None, None> {
    public:
        typedef utopia::None Type;
    };

    template<class Any>
    class FindFunctionSpace {
    public:
        typedef utopia::None Type;
    };

    template<class Space>
    class FindFunctionSpace< TrialFunction<Space> > {
    public:
        typedef Space Type;
    };

    template<class Space>
    class FindFunctionSpace< TestFunction<Space> > {
    public:
        typedef Space Type;
    };

    template<class Space>
    class FindFunctionSpace< TrialFunction<ProductFunctionSpace<Space> > > {
    public:
        typedef Space Type;
    };

    template<class Space>
    class FindFunctionSpace< TestFunction<ProductFunctionSpace<Space> > > {
    public:
        typedef Space Type;
    };

    template<class Left, class Right, class Op>
    class FindFunctionSpace< Binary<Left, Right, Op> > {
    public:
        typedef typename SelectFunctionSpace<
                        typename FindFunctionSpace<Left>::Type,
                        typename FindFunctionSpace<Right>::Type
                        >::Type Type;
    };

    template<class Left, class Right>
    class FindFunctionSpace< Multiply<Left, Right> > {
    public:
        typedef typename SelectFunctionSpace<
                        typename FindFunctionSpace<Left>::Type,
                        typename FindFunctionSpace<Right>::Type
                        >::Type Type;
    };

    template<class Left, class Right>
    class FindFunctionSpace< Equality<Left, Right> > {
    public:
        typedef typename SelectFunctionSpace<
                        typename FindFunctionSpace<Left>::Type,
                        typename FindFunctionSpace<Right>::Type
                        >::Type Type;
    };

    template<class Expr, class Op>
    class FindFunctionSpace< Unary<Expr, Op> > {
    public:
        typedef typename FindFunctionSpace<Expr>::Type Type;
    };

    template<class Expr>
    class FindFunctionSpace< Negate<Expr> > {
    public:
        typedef typename FindFunctionSpace<Expr>::Type Type;
    };

    template<class Expr>
    class FindFunctionSpace< Integral<Expr> > {
    public:
        typedef typename FindFunctionSpace<Expr>::Type Type;
    };

    template<class Expr>
    class FindFunctionSpace< Gradient<Expr> > {
    public:
        typedef typename FindFunctionSpace<Expr>::Type Type;
    };

    template<class Expr>
    class FindFunctionSpace< Divergence<Expr> > {
    public:
        typedef typename FindFunctionSpace<Expr>::Type Type;
    };

    template<class Expr>
    class FindFunctionSpace< Curl<Expr> > {
    public:
        typedef typename FindFunctionSpace<Expr>::Type Type;
    };

    template<class Expr>
    class FindFunctionSpace< Transposed<Expr> > {
    public:
        typedef typename FindFunctionSpace<Expr>::Type Type;
    };

    template<class Expr, class Op>
    class FindFunctionSpace< Reduce<Expr, Op> > {
    public:
        typedef typename FindFunctionSpace<Expr>::Type Type;
    };

    template<class Coeff, class Fun>
    class FindFunctionSpace< Interpolate<Coeff, Fun> > {
    public:
        typedef typename FindFunctionSpace<Fun>::Type Type;
    };

    template<class T>
    class FindFunctionSpace< BlockVar<T> > {
    public:
        typedef utopia::None Type;
    };

    template<class Derived, class AssemblyContext>
    inline int functional_order(const Expression<Derived> &expr, const AssemblyContext &ctx)
    {
        return FunctionalTraits<Derived, AssemblyContext>::order(expr.derived(), ctx);
    }

    template<class Derived, class AssemblyContext>
    inline int functional_type(const Expression<Derived> &expr, const AssemblyContext &ctx)
    {
        return FunctionalTraits<Derived, AssemblyContext>::type(expr.derived(), ctx);
    }
}

#endif //UTOPIA_FE_FUNCTIONAL_TRAITS_HPP
