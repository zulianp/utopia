//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_FACTORY_HPP
#define UTOPIA_UTOPIA_EVAL_FACTORY_HPP

#include "utopia_Eval_Empty.hpp"
// #include "utopia_Tracer.hpp"

namespace utopia {

    // template<class Left, class Right, int Order, class Traits, int Backend>
    // class Eval< Assign<View<Left>, Factory<Right, Order> >, Traits, Backend> {
    // public:
    //     inline static void apply(const Assign<View<Left>, Factory<Right, Order> > &expr)
    //     {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         const auto &left = expr.left();
    //         UTOPIA_BACKEND(Traits).assign_to_range(
    //                 Eval<Left, Traits>::apply(left.expr()),
    //                 expr.right().type(),
    //                 row_range(left),
    //                 col_range(left)
    //         );


    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    // template<class Left, class Right, int Order, class Traits, int Backend>
    // class Eval< Construct<View<Left>, Factory<Right, Order> >, Traits, Backend> {
    // public:
    //     inline static void apply(const Construct<View<Left>, Factory<Right, Order> > &expr)
    //     {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         const auto &left = expr.left();
    //         UTOPIA_BACKEND(Traits).assign_to_range(
    //                 Eval<Left, Traits>::apply(left.expr()),
    //                 expr.right().type(),
    //                 row_range(left),
    //                 col_range(left)
    //         );

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //NEW
    template<class Left, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Factory<LocalZeros, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<LocalZeros, Order>  > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).local_zeros( expr.right().size() );

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Factory<LocalValues<Right>, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<LocalValues<Right>, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).local_values(
                expr.right().size(),
                expr.right().type().value()
            );

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Factory<Values<Right>, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<Values<Right>, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            // UTOPIA_BACKEND(Traits).build(
            //         Eval<Left, Traits>::apply(expr.left()),
            //         expr.right().size(),
            //         expr.right().type()
            // );

            Eval<Left, Traits>::apply(expr.left()).values(
                expr.right().size(),
                expr.right().type().value()
            );

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Factory<NNZ<Right>, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<NNZ<Right>, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).sparse(
                expr.right().size(),
                expr.right().type().nnz()
            );

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Factory<LocalNNZ<Right>, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<LocalNNZ<Right>, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).local_sparse(
                expr.right().size(),
                expr.right().type().nnz()
            );

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Factory<Identity, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<Identity, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).identity(expr.right().size());

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Binary<Number<Right>, Factory<Identity, Order>, Multiplies>>, Traits, Backend> {
    public:
        using Expr = utopia::Assign<Left, Binary<Number<Right>, Factory<Identity, Order>, Multiplies>>;

        inline static void apply(const Expr &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).identity(expr.right().right().size(), expr.right().left());

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Factory<Zeros, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<Zeros, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).zeros(expr.right().size());

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, int Order, class Traits, int Backend>
    class Eval<Construct<Left, Factory<LocalZeros, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Construct<Left, Factory<LocalZeros, Order>  > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).local_zeros( expr.right().size() );

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Construct<Left, Factory<LocalValues<Right>, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Construct<Left, Factory<LocalValues<Right>, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).local_values(
                expr.right().size(),
                expr.right().type().value()
            );

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Construct<Left, Factory<Values<Right>, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Construct<Left, Factory<Values<Right>, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).values(
                expr.right().size(),
                expr.right().type().value()
            );

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, int Order, class Traits, int Backend>
    class Eval<Construct<Left, Factory<Identity, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Construct<Left, Factory<Identity, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).identity(expr.right().size());

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, int Order, class Traits, int Backend>
    class Eval<Construct<Left, Factory<DenseIdentity, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Construct<Left, Factory<DenseIdentity, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).dense_identity(expr.right().size());

            UTOPIA_TRACE_END(expr);
        }
    };


    //NEW
    template<class Left, int Order, class Traits, int Backend>
    class Eval<Assign<Left, Factory<DenseIdentity, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Assign<Left, Factory<DenseIdentity, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).dense_identity(expr.right().size());

            UTOPIA_TRACE_END(expr);
        }
    };

    //NEW
    template<class Left, int Order, class Traits, int Backend>
    class Eval<Construct<Left, Factory<Zeros, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Construct<Left, Factory<Zeros, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            Eval<Left, Traits>::apply(expr.left()).zeros(expr.right().size());

            UTOPIA_TRACE_END(expr);
        }
    };

   
   //NEW
   template<class Left, class Right, int Order, class Traits, int Backend>
   class Eval<Construct<Left, Factory<NNZ<Right>, Order> >, Traits, Backend> {
   public:
       inline static void apply(const Construct<Left, Factory<NNZ<Right>, Order> > &expr)
       {
           UTOPIA_TRACE_BEGIN(expr);

           Eval<Left, Traits>::apply(expr.left()).sparse(
               expr.right().size(),
               expr.right().type().nnz()
           );

           UTOPIA_TRACE_END(expr);
       }
   };

   //NEW
   template<class Left, class Right, int Order, class Traits, int Backend>
   class Eval<Construct<Left, Factory<LocalNNZ<Right>, Order> >, Traits, Backend> {
   public:
       inline static void apply(const Construct<Left, Factory<LocalNNZ<Right>, Order> > &expr)
       {
           UTOPIA_TRACE_BEGIN(expr);

           Eval<Left, Traits>::apply(expr.left()).local_sparse(
               expr.right().size(),
               expr.right().type().nnz()
           );

           UTOPIA_TRACE_END(expr);
       }
   };

   //NEW
   template<class Left, class Right, int Order, class Traits, int Backend>
   class Eval<Construct<Left, Binary<Number<Right>, Factory<Identity, Order>, Multiplies>>, Traits, Backend> {
   public:
       using Expr = utopia::Construct<Left, Binary<Number<Right>, Factory<Identity, Order>, Multiplies>>;

       inline static void apply(const Expr &expr)
       {
           UTOPIA_TRACE_BEGIN(expr);

           Eval<Left, Traits>::apply(expr.left()).identity(expr.right().right().size(), expr.right().left());

           UTOPIA_TRACE_END(expr);
       }
   };

   //NEW
   template<class Left, class Right, int Order, class Traits, int Backend>
   class Eval<Construct<Left, Binary<Number<Right>, Factory<DenseIdentity, Order>, Multiplies>>, Traits, Backend> {
   public:
       using Expr = utopia::Construct<Left, Binary<Number<Right>, Factory<DenseIdentity, Order>, Multiplies>>;

       inline static void apply(const Expr &expr)
       {
           UTOPIA_TRACE_BEGIN(expr);

           Eval<Left, Traits>::apply(expr.left()).dense_identity(expr.right().right().size(), expr.right().left());

           UTOPIA_TRACE_END(expr);
       }
   };


   //NEW
   template<class Expr, class Traits, int Backend>
   class Eval<Factory<Values<Expr>, 1>, Traits, Backend> {
   public:
        using Result = typename Traits::Vector;

       inline static Result apply(const Factory<Values<Expr>, 1> &expr)
       {
           Result result;
           UTOPIA_TRACE_BEGIN(expr);

           result.values(
               expr.size(),
               expr.type().value()
           );

           UTOPIA_TRACE_END(expr);

           return result;
       }
   };

   //NEW
   template<class Expr, int Order, class Traits, int Backend>
   class Eval<Factory<LocalValues<Expr>, Order>, Traits, Backend> {
   public:
        using Result = typename Traits::Vector;

       inline static Result apply(const Factory<LocalValues<Expr>, Order> &expr)
       {
           Result result;
           UTOPIA_TRACE_BEGIN(expr);

           result.local_values(
               expr.size(),
               expr.type().value()
           );

           UTOPIA_TRACE_END(expr);

           return result;
       }
   };

   //NEW
   template<class Traits, int Backend>
   class Eval<Factory<LocalIdentity, 2>, Traits, Backend> {
   public:
        using Result = typename Traits::Matrix;

       inline static Result apply(const Factory<LocalIdentity, 2> &expr)
       {
           Result result;
           // UTOPIA_TRACE_BEGIN(expr);

           result.local_identity(
               expr.size()
           );

           // UTOPIA_TRACE_END(expr);

           return result;
       }
   };

   //NEW FIXME this is for the expression dense()
   template<class Left, int Order, class Traits, int Backend>
   class Eval<Assign<Left, Factory<Resize, Order> >, Traits, Backend> {
   public:
       inline static void apply(const Assign<Left, Factory<Resize, Order> > &expr)
       {
           UTOPIA_TRACE_BEGIN(expr);

           Eval<Left, Traits>::apply(expr.left()).zeros(
               expr.right().size()
           );

           UTOPIA_TRACE_END(expr);
       }
   };

   template<class Left, int Order, class Traits, int Backend>
   class Eval<Construct<Left, Factory<Resize, Order> >, Traits, Backend> {
   public:
       inline static void apply(const Construct<Left, Factory<Resize, Order> > &expr)
       {
           UTOPIA_TRACE_BEGIN(expr);

           Eval<Left, Traits>::apply(expr.left()).zeros(
               expr.right().size()
           );

           UTOPIA_TRACE_END(expr);
       }
   };

}

#endif //UTOPIA_UTOPIA_EVAL_FACTORY_HPP
