#ifndef UTOPIA_UTOPIA_FORWARD_DECLARATIONS_HPP
#define UTOPIA_UTOPIA_FORWARD_DECLARATIONS_HPP

#include <string>

namespace utopia {
    template <class Derived>
    class Expression;

    template <class Derived, int Order>
    class Tensor;

    template <typename T, int BackendType>
    class Backend;

    template <class Left, class Right, class Operation>
    class Binary;

    template <class Expr, class Operation>
    class Unary;

    template <class Expr>
    class Negate;

    template <class Expr>
    class Inverse;

    template <class Expr>
    class Transposed;

    template <typename T>
    class Number;

    template <class Left, class Right>
    class Multiply;

    template <class T>
    class Traits;

    template <class Expr, int Order>
    class Norm;

    template <class Result, int BAKEND_FLAG>
    class Evaluator;

    template <class Left, class Right, class Operation>
    class InPlace;

    template <class Expr>
    class View;

    template <class Expr, int Order>
    class Select;

    template <class Expr>
    class Diag;

    template <class Expr>
    class LocalDiagBlock;

    template <class Expr>
    std::string GetClass();

    template <class Expr>
    class Write;

    template <class Expr>
    class Read;

    template <class Expr, class Operation>
    class Reduce;

    template <class Expr>
    class Trace;

    template <class Left, class Right>
    class Assign;

    // template<class Left, class Right>
    // class Construct;

    template <class Left, class Right>
    using Construct = utopia::Assign<Left, Right>;

    template <class Left, class Right>
    class Equality;

    template <class Left, class Right>
    Construct<Left, Right> construct(Expression<Left> &, const Expression<Right> &);

    class Range;

    template <typename Expr, int Order>
    class Evaluate;

    template <class Expr>
    class Differentiable;

    template <class Expr>
    constexpr int is_differentiable();

    template <class Expr>
    class TreeProperties;

    template <typename T, class Derived>
    T scalar_cast(const Expression<Derived> &);

    template <class Left, class Right>
    class LocalRedistribute;

    template <class Left, class Right>
    class MatrixPtAPProduct;

    template <class Expr, int Number>
    class Variable;

    template <class Expr, class Traits, int Backend>
    class Eval;

    template <class Epxr, class Operation>
    class TensorReduce;

    template <class Type, int Order>
    class Factory;

    template <class Index>
    class Ghosts;

    template <class Type, int Order>
    class SymbolicTensor;

    class Resize;

    template <class Expr>
    class Determinant;

#ifdef WITH_OPENCL
    template <class Expr, int Order>
    class Evaluate;
#endif  // WITH_OPENCL

    template <class Matrix, class Vector, int Backend>
    class QuadraticFunction;

    class SymbolicFunction;

    class Input;

    template <class Vector, int Backend>
    class EvalDots;

    template <class Vector, int Backend>
    class EvalNorm2s;

    template <class Left, class Right, class Traits, int Backend>
    class EvalAssignToView;

    template <typename Scalar, typename SizeType>
    class DistributedMatrix;

    template <typename Scalar, typename SizeType>
    class DistributedVector;

    template <class T>
    class BLAS1Tensor;

    template <class M, class V>
    class BLAS2Matrix;

    template <class M>
    class BLAS3Matrix;

    template <typename Scalar, typename SizeType>
    class AbstractVector;

    template <typename Scalar, typename SizeType>
    class AbstractMatrix;

    template <class V>
    class Operator;

    class Communicator;
    class SelfCommunicator;

    template <class T>
    class RangeDevice;
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_FORWARD_DECLARATIONS_HPP
