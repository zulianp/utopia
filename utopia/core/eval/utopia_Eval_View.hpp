#ifndef UTOPIA_EVAL_VIEW_HPP
#define UTOPIA_EVAL_VIEW_HPP

namespace utopia {

//TODO
// template<class Left, class Right, class Traits, int Backend>
// class Eval< Assign< View<Left>, Right>, Traits, Backend> {
// public:
//     inline static bool apply(const Assign<View<Left>, Right> &expr)
//     {
//         UTOPIA_TRACE_BEGIN(expr);

//         const auto &left = expr.left();
//         auto rr = row_range(left);
//         auto cr = col_range(left);

//         UTOPIA_BACKEND(Traits).assign_to_range(Eval<Left,  Traits>::apply(expr.left().expr()),
//                                                Eval<Right, Traits>::apply(expr.right()),
//                                                rr, cr);
//         UTOPIA_TRACE_END(expr);
//         return true;
//     }
// };

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


//TODO
// template<class Left, class Right, class Traits, int Backend>
// class Eval< Assign<Left, View<Right> >, Traits, Backend> {
// public:
//     inline static bool apply(const Assign<Left, View<Right> > &expr)
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
