#ifndef UTOPIA_EVAL_VIEW_HPP
#define UTOPIA_EVAL_VIEW_HPP

#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

    template<class Left, class Right, class Traits, int Backend>
    class EvalAssignToView {
    public:
        inline static bool apply(const Assign<View<Left>, Right> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            const auto &left = expr.left();
            auto rr = row_range(left);
            auto cr = col_range(left);

            apply_aux(
                Eval<Left,  Traits>::apply(expr.left().expr()),
                rr,
                cr,
                Eval<Right, Traits>::apply(expr.right())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }

        template<class Derived>
        inline static void apply_aux(
            Tensor<Derived, 1> &left,
            const Range &row_range,
            const Range &,
            const Tensor<Derived, 1> &right
        )
        {
            left.derived().assign(row_range, right);
        }

        template<class Derived>
        inline static void apply_aux(
            Tensor<Derived, 2> &left,
            const Range &row_range,
            const Range &col_range,
            const Tensor<Derived, 2> &right
        )
        {
            left.derived().assign(row_range, col_range, right);
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Assign< View<Left>, Right>, Traits, Backend> : public EvalAssignToView<Left, Right, Traits, Backend> {};


//TODO
// template<class Left, class Right, class Traits, int Backend>
// class Eval< Assign< View< Tensor<Left, 1> >, Right>, Traits, Backend> {
// public:
//     typedef utopia::Tensor<Left, 1> LeftTensor;

//     inline static bool apply(const Assign<View<LeftTensor>, Right> &expr)
//     {
//         UTOPIA_TRACE_BEGIN(expr);

//         const auto &left = expr.left();
//         auto rr = row_range(left);
//         auto cr = col_range(left);

//         UTOPIA_BACKEND(Traits).assign_to_range(Eval<LeftTensor, Traits>::apply(expr.left().expr()),
//                                                Eval<Right, Traits>::apply(expr.right()),
//                                                rr, cr);

//         UTOPIA_TRACE_END(expr);
//         return true;
//     }
// };

    template<class Left, int Order, class Right, class Traits, int Backend>
    class Eval< Assign<Tensor<Left, Order>, View<Right> >, Traits, Backend> {
    public:
        inline static bool apply(const Assign<Tensor<Left, Order>, View<Right>> &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            apply_aux(
                Eval<Tensor<Left, Order>, Traits>::apply(expr.left()),
                row_range(expr.right()),
                col_range(expr.right()),
                Eval<Right, Traits>::apply(expr.right().expr())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }

        template<class Derived>
        inline static void apply_aux(
            Tensor<Derived, 1> &left,
            const Range &row_range,
            const Range &,
            const Tensor<Derived, 1> &right
        )
        {
            right.derived().select(row_range, left.derived());
        }

        template<class Derived>
        inline static void apply_aux(
            Tensor<Derived, 2> &left,
            const Range &row_range,
            const Range &col_range,
            const Tensor<Derived, 2> &right
        )
        {
            right.derived().select(row_range, col_range, left.derived());
        }
    };

    template<class Left, class Right, class Traits, int Backend>
    class Eval< Construct<Left, View<Right> >, Traits, Backend> {
    public:
        inline static bool apply(const Construct<Left, View<Right> > &expr)
        {
            UTOPIA_TRACE_BEGIN(expr);

            apply_aux(
                Eval<Left, Traits>::apply(expr.left()),
                row_range(expr.right()),
                col_range(expr.right()),
                Eval<Right, Traits>::apply(expr.right().expr())
            );

            UTOPIA_TRACE_END(expr);
            return true;
        }

        template<class Derived>
        inline static void apply_aux(
            Tensor<Derived, 1> &left,
            const Range &row_range,
            const Range &,
            const Tensor<Derived, 1> &right
        )
        {
            right.derived().select(row_range, left.derived());
        }

        template<class Derived>
        inline static void apply_aux(
            Tensor<Derived, 2> &left,
            const Range &row_range,
            const Range &col_range,
            const Tensor<Derived, 2> &right
        )
        {
            right.derived().select(row_range, col_range, left.derived());
        }
    };


	///////TODOs

	// template<class Left, class Right, class Traits, int Backend>
	// class Eval< Construct< View<Left>, Right>, Traits, Backend> {
	// public:
	//     inline static bool apply(const Construct<View<Left>, Right> &expr)
	//     {
	//         UTOPIA_TRACE_BEGIN(expr);

	//         const auto &left = expr.left();
	//         auto rr = row_range(left);
	//         auto cr = col_range(left);

	//         UTOPIA_BACKEND(Traits).assign_to_range(Eval<Left,  Traits>::apply(expr.left().expr()),
	//                                              Eval<Right, Traits>::apply(expr.right()),
	//                                              rr, cr);

	//         UTOPIA_TRACE_END(expr);
	//         return true;
	//     }
	// };

	// template<class Left, class Right, class Traits, int Backend>
	// class Eval< Construct< View< Tensor<Left, 1> >, Right>, Traits, Backend> {
	// public:
	// typedef utopia::Tensor<Left, 1> LeftTensor;

	//     inline static bool apply(const Construct<View<LeftTensor>, Right> &expr)
	//     {
	//         UTOPIA_TRACE_BEGIN(expr);

	//         const auto &left = expr.left();
	//         auto rr = row_range(left);
	//         auto cr = col_range(left);

	//         UTOPIA_BACKEND(Traits).assign_to_range(Eval<LeftTensor, Traits>::apply(expr.left().expr()),
	//                                              Eval<Right, Traits>::apply(expr.right()),
	//                                              rr, cr);

	//         UTOPIA_TRACE_END(expr);
	//         return true;
	//     }
	// };


	// template<class Left, class Right, class Traits, int Backend>
	// class Eval< Construct<Left, View<Right> >, Traits, Backend> {
	// public:
	//     inline static bool apply(const Construct<Left, View<Right> > &expr)
	//     {
	//         UTOPIA_TRACE_BEGIN(expr);

	//         UTOPIA_BACKEND(Traits).assign_from_range(
	//                 Eval<Left,  Traits>::apply(expr.left()),
	//                 Eval<Right, Traits>::apply(expr.right().expr()),
	//                 row_range(expr.right()),
	//                 col_range(expr.right())
	//         );

	//         UTOPIA_TRACE_END(expr);
	//         return true;
	//     }
	// };

}

#endif //UTOPIA_EVAL_VIEW_HPP
