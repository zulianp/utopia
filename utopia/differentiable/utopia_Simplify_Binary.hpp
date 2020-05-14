#ifndef DO_NOT_SIMPLIFY_DERIVATIVES

#ifndef UTOPIA_SIMPLIFY_BINARY_HPP
#define UTOPIA_SIMPLIFY_BINARY_HPP

namespace utopia {

    template<class Left, class Right, class Operation>
    class BinarySimplify {
    public:
        typedef utopia::Binary<Left, Right, Operation> Type;

        static UTOPIA_STORE_CONST(Type) make(const Left &left, const Right &right, const Operation &op)
        {
            return Type(left, right, op);
        }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////// PLUS ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    template<class Left, int Order>
    class BinarySimplify<Left, Factory<Zeros, Order>, Plus> {
    public:
        using Type = Left;

        static UTOPIA_STORE_CONST(Type) make(const Left &left, const Factory<Zeros, Order> &, const Plus &)
        {
            return left;
        }
    };

    template<int Order>
    class BinarySimplify<Factory<Identity, Order>, Factory<Identity, Order>, Plus> {
    public:
        typedef typename utopia::Factory<Identity, Order>::Scalar Scalar;
        typedef utopia::Binary< Number<Scalar>, Factory<Identity, Order>, Multiplies> Type;

        static UTOPIA_STORE_CONST(Type) make(const Factory<Identity, Order> &left, const Factory<Identity, Order> &, const Plus &)
        {
            return Scalar(2) * left;
        }
    };

    template<int Order>
    class BinarySimplify<Factory<Zeros, Order>, Factory<Zeros, Order>, Plus> {
    public:
        typedef Factory<Zeros, Order> Type;

        static UTOPIA_STORE_CONST(Type) make(const Factory<Zeros, Order> &left, const Factory<Zeros, Order> &, const Plus &)
        {
            return left;
        }
    };

    template< int Order, class Right>
    class BinarySimplify<Factory<Zeros, Order>, Right, Plus> {
    public:
        using Type = Right;

        static UTOPIA_STORE_CONST(Type) make(const Factory<Zeros, Order> &, const Right &right, const Plus &)
        {
            return right;
        }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////// EMULTIPLIES ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    template<class Expr, int Order>
    class BinarySimplify<Expr, Factory<Identity, Order>, EMultiplies> {
    public:
        using Type = typename Simplify<Expr>::Type;

        static const UTOPIA_STORE_CONST(Type) make(const Expr &left, const Factory<Identity, Order> &/*right*/, const EMultiplies &)
        {
            return Simplify<Expr>::make(left);
        }
    };

    template<class Expr, int Order>
    class BinarySimplify< Factory<Identity, Order>, Expr, EMultiplies> {
    public:
        using Type = typename Simplify<Expr>::Type;

        static UTOPIA_STORE_CONST(Type) make(const Factory<Identity, Order> &/*left*/, const Expr &right, const EMultiplies &)
        {
            return Simplify<Expr>::make(right);
        }
    };

    template<class Expr, int Order>
    class BinarySimplify<Factory<Zeros, Order>, Expr, EMultiplies>{
    public:
        typedef utopia::Factory<Zeros, Order> Type;

        static UTOPIA_STORE_CONST(Type) make(const Factory<Zeros, Order> &left, const Expr & /*right*/, const EMultiplies &)
        {
            return left;
        }
    };

    template<class Expr, int Order>
    class BinarySimplify<Expr, Factory<Zeros, Order>, EMultiplies>{
    public:
        typedef utopia::Factory<Zeros, Order> Type;

        static UTOPIA_STORE_CONST(Type) make(const Expr & /*left*/, const Factory<Zeros, Order> &right, const EMultiplies &)
        {
            return right;
        }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////// MULTIPLIES ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    template<class Expr, int Order>
    class BinarySimplify<Expr, Factory<Identity, Order>, Multiplies> {
    public:
        using Type = typename Simplify<Expr>::Type;

        static const UTOPIA_STORE_CONST(Type) make(const Expr &left, const Factory<Identity, Order> &/*right*/, const Multiplies &)
        {
            return Simplify<Expr>::make(left);
        }
    };

    template<class Expr, int Order>
    class BinarySimplify< Factory<Identity, Order>, Expr, Multiplies> {
    public:
        using Type = typename Simplify<Expr>::Type;

        static UTOPIA_STORE_CONST(Type) make(const Factory<Identity, Order> &/*left*/, const Expr &right, const Multiplies &)
        {
            return Simplify<Expr>::make(right);
        }
    };

    template<class Expr, int Order>
    class BinarySimplify<Factory<Zeros, Order>, Expr, Multiplies>{
    public:
        typedef utopia::Factory<Zeros, Order> Type;

        static UTOPIA_STORE_CONST(Type) make(const Factory<Zeros, Order> &left, const Expr & /*right*/, const Multiplies &)
        {
            return left;
        }
    };

    template<class Expr, int Order>
    class BinarySimplify<Expr, Factory<Zeros, Order>, Multiplies>{
    public:
        typedef utopia::Factory<Zeros, Order> Type;

        static UTOPIA_STORE_CONST(Type) make(const Expr & /*left*/, const Factory<Zeros, Order> &right, const Multiplies &)
        {
            return right;
        }
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template<class Left, class Right, class Operation>
    class Simplify< Binary<Left, Right, Operation> >{
    public:
        using SLeft = typename Simplify<Left>::Type;
        using SRight = typename Simplify<Right>::Type;

        typedef typename utopia::BinarySimplify<SLeft, SRight, Operation> Sim;
        using Type = typename Sim::Type;

        static UTOPIA_STORE_CONST(Type) make(const Binary<Left, Right, Operation> &expr)
        {
            return Sim::make(
                             Simplify<Left>::make(expr.left()),
                             Simplify<Right>::make(expr.right()),
                             expr.operation() );
        }
    };
}

#endif //UTOPIA_SIMPLIFY_BINARY_HPP
#endif //DO_NOT_SIMPLIFY_DERIVATIVES
