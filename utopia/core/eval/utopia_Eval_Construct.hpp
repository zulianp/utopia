// #ifndef UTOPIA_UTOPIA_EVAL_CONSTRUCT_HPP_HPP
// #define UTOPIA_UTOPIA_EVAL_CONSTRUCT_HPP_HPP

// #include "utopia_ForwardDeclarations.hpp"
// #include "utopia_Eval_Empty.hpp"
// #include "utopia_Tracer.hpp"

// namespace utopia {
    
//     template<class Left, class RowPtr, class ColIndex, class Values, class Traits, int Backend>
//     class Eval< Construct<
//                         Left,
//                         Factory< CRS<RowPtr, ColIndex, Values>, 2>>,
//                         Traits,
//                         Backend
//                         > {
//     public:
//         inline static bool apply(const Construct<Left, Factory<CRS<RowPtr, ColIndex, Values>, 2>> &expr)
//         {
//             UTOPIA_TRACE_BEGIN(expr);
//             auto && left = Eval<Left, Traits>::apply(expr.left());
//             auto && crs  = expr.right().type();
//             auto s       = expr.right().size();

//             //FIXME
//             left.crs_init(
//                 left.comm().get(),
//                 INVALID_INDEX,
//                 INVALID_INDEX,
//                 s.get(0),
//                 s.get(1),
//                 crs.rowPtr(),
//                 crs.cols(),
//                 crs.values()
//             );

//             UTOPIA_TRACE_END(expr);
//             return true;
//         }
//     };

// }

// #endif //UTOPIA_UTOPIA_EVAL_CONSTRUCT_HPP_HPP
