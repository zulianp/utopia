#ifndef UTOPIA_TRAVERSE_HPP
#define UTOPIA_TRAVERSE_HPP

#include "utopia_FEForwardDeclarations.hpp"
#include <iostream>
#include <memory>

namespace utopia {
    static const int TRAVERSE_CONTINUE = 0;
    static const int TRAVERSE_STOP = 1;
    static const int TRAVERSE_SKIP_SUBTREE = 2;


    template<class Expr, class Visitor>
    inline static int traverse(Expr &expr, Visitor &visitor)
    {
        std::cout << "[Warning] Traverse: encountered unhandled expression: " << expr.getClass() << std::endl;
        return TRAVERSE_CONTINUE;
    }

    template<class Out, class F, class Visitor>
    inline static int traverse(const ContextFunction<Out, F> &expr, Visitor &visitor)
    {
        visitor.visit(expr);
        return TRAVERSE_CONTINUE;
    }

    template<class Visitor>
    inline static int traverse(const SymbolicFunction &expr, Visitor &visitor)
    {
        visitor.visit(expr);
        return TRAVERSE_CONTINUE;
    }


    template<class Type, int Order, class Visitor>
    inline static int traverse(const SymbolicTensor<Type, Order> &expr, Visitor &visitor)
    {
        visitor.visit(expr);
        return TRAVERSE_CONTINUE;
    }


    template<class Type, int Order, class Visitor>
    inline static int traverse(const Factory<Type, Order> &expr, Visitor &visitor)
    {
        visitor.visit(expr);
        return TRAVERSE_CONTINUE;
    }

    template<class Expr, class Visitor>
    inline static int traverse(Integral<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(Trace<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }


    template<class Expr, class Visitor>
    inline static int traverse(Negate<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr)) {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }


    template<class Expr, class Visitor>
    inline static int traverse(Determinant<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr)) {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(Inverse<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr)) {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Coefficient, class Function, class Visitor>
    inline static int traverse(Interpolate<Coefficient, Function> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class Expr, class Operation, class Visitor>
    inline static int traverse(Unary<Expr, Operation> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr)) {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Operation, class Visitor>
    inline static int traverse(Reduce<Expr, Operation> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(Divergence<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(Transposed<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(Gradient<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(Curl<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class T, int Order, class Visitor>
    inline static int traverse(Wrapper<T, Order> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class T, class Visitor>
    inline static int traverse(Number<T> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class T, int Order, class Visitor>
    inline static int traverse(ConstantCoefficient<T, Order> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class Left, class Right, class Operation, class Visitor>
    inline static int traverse(Binary<Left, Right, Operation> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                if(traverse(expr.left(), visitor) == TRAVERSE_CONTINUE) {
                    return traverse(expr.right(), visitor);
                }
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }

    }

    template<class Left, class Right, class Operation, class Visitor>
    inline static int traverse(Binary<Left, Number<Right>, Operation> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.left(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Visitor>
    inline static int traverse(SymbolicFunction &expr, Visitor &visitor)
    {
        visitor.visit(expr);
        return TRAVERSE_CONTINUE;
    }

    template<class Left, class Right, class Operation, class Visitor>
    inline static int traverse(Binary<Left, BlockVar<Right>, Operation> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.left(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Left, class Right, class Visitor>
    inline static int traverse(Multiply<Left, Right> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                if(traverse(expr.left(), visitor) == TRAVERSE_CONTINUE) {
                    return traverse(expr.right(), visitor);
                }
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Left, class Right, class Visitor>
    inline static int traverse(Equality<Left, Right> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                if(traverse(expr.left(), visitor) == TRAVERSE_CONTINUE) {
                    return traverse(expr.right(), visitor);
                }
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class FunctionSpaceT, class Visitor>
    inline static int traverse(TrialFunction<FunctionSpaceT> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class FunctionSpaceT, class Visitor>
    inline static int traverse(TestFunction<FunctionSpaceT> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }


    //const versions

    template<class T, class Visitor>
    inline static int traverse(const BlockVar<T> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }


    template<class Expr, class Visitor>
    inline static int traverse(const Determinant<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr)) {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(const Inverse<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr)) {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(const Expr &expr, Visitor &visitor)
    {
        std::cout << "[Warning] Encountered unhandled expression: " << expr.getClass() << std::endl;
        return TRAVERSE_CONTINUE;
    }

    template<class Expr, class Visitor>
    inline static int traverse(const Integral<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }


    template<class Expr, class Visitor>
    inline static int traverse(const Trace<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(const Negate<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr)) {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }


    template<class Expr, class Visitor>
    inline static int traverse(Dual<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr)) {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(const Dual<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr)) {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Coefficient, class Function, class Visitor>
    inline static int traverse(const Interpolate<Coefficient, Function> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class Expr, class Operation, class Visitor>
    inline static int traverse(const Unary<Expr, Operation> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr)) {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Operation, class Visitor>
    inline static int traverse(const Reduce<Expr, Operation> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(const Divergence<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(const Transposed<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(const Gradient<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Expr, class Visitor>
    inline static int traverse(const Curl<Expr> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.expr(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class T, int Order, class Visitor>
    inline static int traverse(const Wrapper<T, Order> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class T, class Visitor>
    inline static int traverse(const Number<T> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class T, int Order, class Visitor>
    inline static int traverse(const ConstantCoefficient<T, Order> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class Left, class Right, class Operation, class Visitor>
    inline static int traverse(const Binary<Left, Right, Operation> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                if(traverse(expr.left(), visitor) == TRAVERSE_CONTINUE) {
                    return traverse(expr.right(), visitor);
                }
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }

    }

    template<class Left, class Right, class Operation, class Visitor>
    inline static int traverse(const Binary<Left, Number<Right>, Operation> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.left(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Left, class Right, class Operation, class Visitor>
    inline static int traverse(const Binary<Left, BlockVar<Right>, Operation> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.left(), visitor);
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Left, class Right, class Visitor>
    inline static int traverse(const Multiply<Left, Right> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                if(traverse(expr.left(), visitor) == TRAVERSE_CONTINUE) {
                    return traverse(expr.right(), visitor);
                }
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class Left, class Right, class Visitor>
    inline static int traverse(const Equality<Left, Right> &expr, Visitor &visitor)
    {
        switch(visitor.visit(expr))
        {
            case TRAVERSE_CONTINUE:
            {
                if(traverse(expr.left(), visitor) == TRAVERSE_CONTINUE) {
                    return traverse(expr.right(), visitor);
                }
            }

            case TRAVERSE_STOP:
            {
                return TRAVERSE_STOP;
            }

            case TRAVERSE_SKIP_SUBTREE:
            {
                return TRAVERSE_CONTINUE;
            }

            default: {
                std::cout << "[Error] INVALID RETURN VALUE: stopping traversal" << std::endl;
                return TRAVERSE_STOP;
            }
        }
    }

    template<class FunctionSpaceT, class Visitor>
    inline static int traverse(const TrialFunction<FunctionSpaceT> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class FunctionSpaceT, class Visitor>
    inline static int traverse(const TestFunction<FunctionSpaceT> &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class Eq, class Visitor>
    inline static int traverse(const Equations<Eq> &expr, Visitor &visitor)
    {
        const int ret = visitor.visit(expr);
        switch(ret)
        {
            case TRAVERSE_CONTINUE:
            {
                return traverse(expr.template get<0>(), visitor);
            }
            default: {
                return ret;
            }
        }
    }

    template<class Visitor>
    class EquationsTraverse {
    public:
        EquationsTraverse(Visitor &visitor)
        : visitor(visitor), traverse_command(TRAVERSE_CONTINUE) {}

        template<class Eq>
        inline void operator()(const int, Eq &eq) {
            if(traverse_command == TRAVERSE_STOP) return;
            traverse_command = traverse(eq, visitor);
        }

        Visitor &visitor;
        int traverse_command;
    };

    template<class... Eqs, class Visitor>
    inline static int traverse(const Equations<Eqs...> &expr, Visitor &visitor)
    {
        const int ret = visitor.visit(expr);
        switch(ret)
        {
            case TRAVERSE_CONTINUE:
            {
                EquationsTraverse<Visitor> eq_traverse(visitor);
                expr.each(eq_traverse);
                return eq_traverse.traverse_command;
            }

            default:
            {
                return ret;
            }
        }
    }

    template<class Expr>
    class FindExpression {
    public:

        template<class Any>
        inline constexpr static int visit(const Any &) { return TRAVERSE_CONTINUE; }

        int visit(const Expr &expr)
        {
            expr_ = &expr;
            return TRAVERSE_STOP;
        }

        FindExpression()
        : expr_(nullptr)
        {}

        inline bool found() const
        {
            return expr_ != nullptr;
        }

        inline const Expr &get() const
        {
            return *expr_;
        }

        template<class ExprTree>
        inline bool apply(const ExprTree &expr)
        {
            traverse(expr, *this);
            return found();
        }

        const Expr * expr_;
    };


    template<template<class...> class Expr>
    class TPLExpressionExists {
    public:

        template<class Any>
        inline constexpr static int visit(const Any &) { return TRAVERSE_CONTINUE; }

        template<class... Inner>
        int visit(const Expr<Inner...> &expr)
        {
            found_ = true;
            return TRAVERSE_STOP;
        }

        TPLExpressionExists()
        : found_(false)
        {}

        inline bool found() const
        {
            return found_;
        }


        template<class ExprTree>
        inline bool apply(const ExprTree &expr)
        {
            traverse(expr, *this);
            return found();
        }

        bool found_;
    };

    template<class FunctionSpaceT, class ExprTree>
    inline std::shared_ptr<FunctionSpaceT> trial_space(const ExprTree &tree)
    {
        std::shared_ptr<FunctionSpaceT> ret = nullptr;

        FindExpression< TrialFunction<FunctionSpaceT> > f;
        if(f.apply(tree)) {
            ret = f.get().space_ptr();
        }

        return ret;
    }

    template<class FunctionSpaceT, class ExprTree>
    inline std::shared_ptr<FunctionSpaceT> test_space(const ExprTree &tree)
    {
        std::shared_ptr<FunctionSpaceT> ret = nullptr;

        FindExpression< TestFunction<FunctionSpaceT> > f;
        if(f.apply(tree)) {
            ret = f.get().space_ptr();
        }

        return ret;
    }

    template<class ExprTree>
    inline bool is_trial(const ExprTree &tree)
    {
        TPLExpressionExists<utopia::TrialFunction> f;
        return f.apply(tree);
    }

    template<class ExprTree>
    inline bool is_test(const ExprTree &tree)
    {
        TPLExpressionExists<utopia::TestFunction> f;
        return f.apply(tree);
    }
}

#endif //UTOPIA_TRAVERSE_HPP
