//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_FACTORY_HPP
#define UTOPIA_UTOPIA_EVAL_FACTORY_HPP

#include "utopia_Eval_Empty.hpp"

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

            // UTOPIA_BACKEND(Traits).build(
            //         Eval<Left, Traits>::apply(expr.left()),
            //         expr.right().size(),
            //         expr.right().type()
            // );

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



    // template<class Left, class Right, int Order, class Traits, int Backend>
    // class Eval<Assign<Left, Factory<Right, Order> >, Traits, Backend> {
    // public:
    //     inline static void apply(const Assign<Left, Factory<Right, Order> > &expr)
    //     {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         UTOPIA_BACKEND(Traits).build(
    //                 Eval<Left, Traits>::apply(expr.left()),
    //                 expr.right().size(),
    //                 expr.right().type()
    //         );


    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    // template<class Left, class Right, int Order, class Traits, int Backend>
    // class Eval<Construct<Left, Factory<Right, Order> >, Traits, Backend> {
    // public:
    //     inline static void apply(const Construct<Left, Factory<Right, Order> > &expr)
    //     {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         UTOPIA_BACKEND(Traits).build(
    //                 Eval<Left, Traits>::apply(expr.left()),
    //                 expr.right().size(),
    //                 expr.right().type()
    //         );


    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    //NEW
    template<class Left, class Right, int Order, class Traits, int Backend>
    class Eval<Construct<Left, Factory<LocalValues<Right>, Order> >, Traits, Backend> {
    public:
        inline static void apply(const Construct<Left, Factory<LocalValues<Right>, Order> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            // UTOPIA_BACKEND(Traits).build(
            //         Eval<Left, Traits>::apply(expr.left()),
            //         expr.right().size(),
            //         expr.right().type()
            // );

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
}

#endif //UTOPIA_UTOPIA_EVAL_FACTORY_HPP
