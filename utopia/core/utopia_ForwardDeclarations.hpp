//
// Created by Patrick Zulian on 18/05/15.
//

#ifndef utopia_utopia_FORWARDDECLARATIONS_HPP
#define utopia_utopia_FORWARDDECLARATIONS_HPP

#include <string>

namespace utopia {
    template<class Derived>
    class Expression;


    template<typename T, int BackendType>
    class Backend;

    template<class Left, class Right, class Operation>
    class Binary;

    template<class Expr, class Operation>
    class Unary;

    template<typename T>
    class Number;

    template<class Left, class Right>
    class Multiply;

    template<class Implementation, int Order>
    class Wrapper;

    template<class T>
    class Traits;

    template <class T>
    class Matrix;

    template<class Expr, int Order>
    class Norm;

    template <class Result, int BAKEND_FLAG>
    class Evaluator;

    template<class Left, class Right, class Operation>
    class InPlace;

    template<class Expr>
    class View;

    template<class Expr, typename SizeType, int Order>
    class Select;

    // template<class Expr>
    // class Diag;

    template<class Expr>
    class LocalDiagBlock;


    template<class Expr>
    std::string GetClass();


    template<class Expr>
    class Write;


    template<class Expr>
    class Read;

    template<class Expr, class Operation>
    class Reduce;

    template<class Expr>
    class Trace;

    template<class Left, class Right>
    class Assign;

    template<class Left, class Right>
    class Construct;

    template<class Left, class Right>
    Construct<Left, Right> construct(Expression<Left> &, const Expression<Right> &);

    class Range;

    template<typename Expr, int Order>
    class Evaluate;

    template<class Expr>
    class Differentiable;

    // template<class Expr>
    // class Derivative;

    template<class Expr>
    constexpr int is_differentiable();


    template<class Expr>
    class TreeProperties;

    template<typename T, class Derived>
    T scalar_cast(const Expression<Derived> &);


    template<class Left, class Right>
    class LocalRedistribute;

    template<class Left, class Right>
    class MatrixPtAPProduct;

    template<class Expr, int Number>
    class Variable;

    template<class Expr, class Traits, int Backend>
    class Eval;

    template<class Epxr, class Operation>
    class TensorReduce;

    template<class Type, int Order>
    class Factory;

    class Resize;

    template<class Expr>
    class Determinant;

#ifdef WITH_OPENCL
    template<class Expr, int Order>
    class Evaluate;
#endif //WITH_OPENCL    

}

#endif //utopia_utopia_FORWARDDECLARATIONS_HPP
