#ifndef UTOPIA_UTOPIA_EVAL_FACTORY_HPP
#define UTOPIA_UTOPIA_EVAL_FACTORY_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    ///////////////////////////////////////////////////////////////////////////////////////////////////

    template <class Left, class Traits, int Backend>
    class Eval<Assign<Left, Factory<LocalZeros, 1>>, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<LocalZeros, 1>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            if (Backend == BLAS) {
                Eval<Left, Traits>::apply(expr.left()).zeros(serial_layout(expr.right().size().get(0)));
            } else {
                Eval<Left, Traits>::apply(expr.left())
                    .zeros(
                        layout(Traits::Communicator::get_default(), expr.right().size().get(0), Traits::determine()));
            }

            UTOPIA_TRACE_END(expr);
        }
    };  // namespace utopia

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Left, Factory<LocalValues<Right>, 1>>, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<LocalValues<Right>, 1>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            if (Backend == BLAS) {
                Eval<Left, Traits>::apply(expr.left())
                    .values(serial_layout(expr.right().size().get(0)), expr.right().type().value());

            } else {
                Eval<Left, Traits>::apply(expr.left())
                    .values(
                        layout(Traits::Communicator::get_default(), expr.right().size().get(0), Traits::determine()),
                        expr.right().type().value());
            }

            UTOPIA_TRACE_END(expr);
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////

    template <class Left, class Traits, int Backend>
    class Eval<Assign<Left, Factory<Zeros, 1>>, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<Zeros, 1>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left())
                .zeros(layout(Traits::Communicator::get_default(), Traits::decide(), expr.right().size().get(0)));

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class Eval<Assign<Left, Factory<Values<Right>, 1>>, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<Values<Right>, 1>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            // Eval<Left, Traits>::apply(expr.left()).values(expr.right().size(), expr.right().type().value());

            Eval<Left, Traits>::apply(expr.left())
                .values(layout(Traits::Communicator::get_default(), Traits::decide(), expr.right().size().get(0)),
                        expr.right().type().value());

            UTOPIA_TRACE_END(expr);
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////

    template <class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Factory<NNZ<Right>, Order>>, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<NNZ<Right>, Order>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            // Eval<Left, Traits>::apply(expr.left()).sparse(expr.right().size(), expr.right().type().nnz());
            Eval<Left, Traits>::apply(expr.left())
                .sparse(layout(Traits::Communicator::get_default(),
                               Traits::decide(),
                               Traits::decide(),
                               expr.right().size().get(0),
                               expr.right().size().get(1)),
                        expr.right().type().nnz(),
                        expr.right().type().nnz());

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Factory<LocalNNZ<Right>, Order>>, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<LocalNNZ<Right>, Order>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            // Eval<Left, Traits>::apply(expr.left()).local_sparse(expr.right().size(), expr.right().type().nnz());
            Eval<Left, Traits>::apply(expr.left())
                .sparse(layout(Traits::Communicator::get_default(),
                               expr.right().size().get(0),
                               expr.right().size().get(1),
                               Traits::determine(),
                               Traits::determine()),
                        expr.right().type().nnz(),
                        expr.right().type().nnz());

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class Left, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Factory<Identity, Order>>, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<Identity, Order>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            // Eval<Left, Traits>::apply(expr.left()).identity(expr.right().size());
            Eval<Left, Traits>::apply(expr.left())
                .identity(layout(Traits::Communicator::get_default(),
                                 Traits::decide(),
                                 Traits::decide(),
                                 expr.right().size().get(0),
                                 expr.right().size().get(1)),
                          1.0);

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class Left, int Order, class Traits, int Backend>
    class Eval<Assign<Tensor<Left, Order>, Factory<DenseIdentity, Order>>, Traits, Backend> {
    public:
        inline static void apply(const Assign<Tensor<Left, Order>, Factory<DenseIdentity, Order>> &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Tensor<Left, Order>, Traits>::apply(expr.left())
                .dense_identity(layout(Traits::Communicator::get_default(),
                                       Traits::decide(),
                                       Traits::decide(),
                                       expr.right().size().get(0),
                                       expr.right().size().get(1)),
                                1.0);

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class Traits, int Backend>
    class Eval<Factory<LocalIdentity, 2>, Traits, Backend> {
    public:
        using Expr = utopia::Factory<LocalIdentity, 2>;
        using Result = typename Traits::Matrix;

        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            result.identity(layout(Traits::Communicator::get_default(),
                                   expr.size().get(0),
                                   expr.size().get(1),
                                   Traits::determine(),
                                   Traits::determine()),
                            1.0);

            UTOPIA_TRACE_END(expr);
        }
    };

    template <class Traits, int Backend>
    class Eval<Factory<LocalDenseIdentity, 2>, Traits, Backend> {
    public:
        using Expr = utopia::Factory<LocalDenseIdentity, 2>;
        using Result = typename Traits::Matrix;
        UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

        inline static void apply(const Expr &expr, Result &result) {
            UTOPIA_TRACE_BEGIN(expr);

            // result.local_dense_identity(expr.size());

            result.identity(layout(Traits::Communicator::get_default(),
                                   Traits::decide(),
                                   Traits::decide(),
                                   expr.size().get(0),
                                   expr.size().get(1)),
                            1.0);

            UTOPIA_TRACE_END(expr);
        }
    };

    //////////////////////////////////////////////////////////////////////////////////////////
    // composites
    //////////////////////////////////////////////////////////////////////////////////////////

    // template <class Left, class Right, int Order, class Traits, int Backend>
    // class Eval<Assign<Tensor<Left, Order>, Binary<Number<Right>, Factory<Identity, Order>, Multiplies>>,
    //            Traits,
    //            Backend> {
    // public:
    //     using Expr = utopia::Assign<Tensor<Left, Order>, Binary<Number<Right>, Factory<Identity, Order>,
    //     Multiplies>>;

    //     inline static void apply(const Expr &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         Eval<Tensor<Left, Order>, Traits>::apply(expr.left())
    //             .identity(expr.right().right().size(), expr.right().left());

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    // template <class Left, class Right, int Order, class Traits, int Backend>
    // class Eval<Construct<Left, Binary<Number<Right>, Factory<Identity, Order>, Multiplies>>, Traits, Backend> {
    // public:
    //     using Expr = utopia::Construct<Left, Binary<Number<Right>, Factory<Identity, Order>, Multiplies>>;

    //     inline static void apply(const Expr &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         Eval<Left, Traits>::apply(expr.left()).identity(expr.right().right().size(), expr.right().left());

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    // template <class Left, class Right, int Order, class Traits, int Backend>
    // class Eval<Construct<Left, Binary<Number<Right>, Factory<LocalIdentity, Order>, Multiplies>>, Traits, Backend> {
    // public:
    //     using Expr = utopia::Construct<Left, Binary<Number<Right>, Factory<LocalIdentity, Order>, Multiplies>>;

    //     inline static void apply(const Expr &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         Eval<Left, Traits>::apply(expr.left()).local_identity(expr.right().right().size(), expr.right().left());

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    // template <class Left, class Right, int Order, class Traits, int Backend>
    // class Eval<Construct<Left, Binary<Number<Right>, Factory<DenseIdentity, Order>, Multiplies>>, Traits, Backend> {
    // public:
    //     using Expr = utopia::Construct<Left, Binary<Number<Right>, Factory<DenseIdentity, Order>, Multiplies>>;

    //     inline static void apply(const Expr &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         Eval<Left, Traits>::apply(expr.left()).dense_identity(expr.right().right().size(), expr.right().left());

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    // template <class InnerExpr, class Traits, int Backend>
    // class Eval<Factory<Values<InnerExpr>, 1>, Traits, Backend> {
    // public:
    //     using Expr = utopia::Factory<Values<InnerExpr>, 1>;
    //     using Result = typename Traits::Vector;

    //     UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

    //     inline static void apply(const Expr &expr, Result &result) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         result.values(expr.size(), expr.type().value());

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    // template <class InnerExpr, int Order, class Traits, int Backend>
    // class Eval<Factory<LocalValues<InnerExpr>, Order>, Traits, Backend> {
    // public:
    //     using Expr = utopia::Factory<LocalValues<InnerExpr>, Order>;
    //     using Result = typename Traits::Vector;

    //     UTOPIA_EVAL_APPLY_TO_TEMPORARY(Expr, Result)

    //     inline static void apply(const Expr &expr, Result &result) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         result.local_values(expr.size(), expr.type().value());

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    // NEW FIXME this is for the expression dense()
    // template <class Left, int Order, class Traits, int Backend>
    // class Eval<Assign<Left, Factory<Resize, Order>>, Traits, Backend> {
    // public:
    //     inline static void apply(const Assign<Left, Factory<Resize, Order>> &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         Eval<Left, Traits>::apply(expr.left()).zeros(expr.right().size());

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    // template <class Derived, typename T>
    // using AssignMatrixWithIdShit =
    //     utopia::Assign<Tensor<Derived, 2>,
    //                    Binary<Tensor<Derived, 2>, Binary<Number<T>, Factory<LocalIdentity, 2>, Multiplies>, Plus>>;

    // template <class Derived, typename T, class Traits, int Backend>
    // class Eval<AssignMatrixWithIdShit<Derived, T>, Traits, Backend> {
    // public:
    //     using Matrix = utopia::Tensor<Derived, 2>;
    //     using Scalar = typename Traits::Scalar;

    //     inline static void apply(const AssignMatrixWithIdShit<Derived, T> &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         auto &&l = Eval<Matrix, Traits>::apply(expr.left());
    //         const Scalar diag = expr.right().right().left();

    //         l.construct(Eval<Matrix, Traits>::apply(expr.right().left()));
    //         l.shift_diag(diag);

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    // template <class Left, class Traits, int Backend>
    // class Eval<InPlace<Left, Factory<Identity, 2>, Minus>, Traits, Backend> {
    // public:
    //     inline static bool apply(const InPlace<Left, Factory<Identity, 2>, Minus> &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         Eval<Left, Traits>::apply(expr.left()).shift_diag(-1.);

    //         UTOPIA_TRACE_END(expr);
    //         return true;
    //     }
    // };

    // template <class Left, class Traits, int Backend>
    // class Eval<InPlace<Left, Factory<Identity, 2>, Plus>, Traits, Backend> {
    // public:
    //     inline static bool apply(const InPlace<Left, Factory<Identity, 2>, Plus> &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         Eval<Left, Traits>::apply(expr.left()).shift_diag(1.);

    //         UTOPIA_TRACE_END(expr);
    //         return true;
    //     }
    // };

}  // namespace utopia

#endif  // UTOPIA_UTOPIA_EVAL_FACTORY_HPP
