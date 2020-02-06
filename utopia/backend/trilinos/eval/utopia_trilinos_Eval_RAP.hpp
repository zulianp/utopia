#ifndef UTOPIA_TRILINOS_EVAL_RAP_HPP
#define UTOPIA_TRILINOS_EVAL_RAP_HPP

#include "utopia_trilinos_ForwardDeclarations.hpp"

//FIXME find right macro
#ifdef WITH_TRILINOS_TPETRAEXT

#include "TpetraExt_TripleMatrixMultiply.hpp"

//useful links:
//https://trilinos.org/docs/dev/packages/tpetra/doc/html/namespaceTpetra_1_1TripleMatrixMultiply.html
//see MultiplyRAP

namespace utopia {

    // mat-mat-mat multiplication
    template<class M1, class M2, class M3, class Traits>
    class Eval< Multiply< Multiply< Transposed<Tensor<M1, 2>>, Tensor<M2, 2> >, Tensor<M3, 2> >, Traits, TRILINOS> {
    public:
        typedef utopia::Multiply< Multiply< Transposed<Tensor<M1, 2>>, Tensor<M2, 2> >, Tensor<M3, 2> > Expr;
        typedef EXPR_TYPE(Traits, Expr) Result;

        inline static Result apply(const Expr &expr)
        {
            Result result;
            apply(expr, result);
            return result;
        }

        inline static void apply(const Expr &expr, Result &result)
        {
            UTOPIA_TRACE_BEGIN(expr);

            try {

                auto &&R = expr.left().left().expr();
                auto &&A = expr.left().right();
                auto &&P = expr.right();

                UTOPIA_REPORT_ALLOC("Eval_RAP");
                auto raw_result = Teuchos::rcp(new typename Result::CrsMatrixType(
                    raw_type(R)->getDomainMap(),
                    0,
                    Tpetra::StaticProfile
                ));

                assert(!empty(R));
                assert(!empty(A));
                assert(!empty(P));

                //Performs optimal triple product
                //Ac = R*A*P,
                Tpetra::TripleMatrixMultiply::MultiplyRAP(
                    *raw_type(R),
                    true, //transposeR
                    *raw_type(A),
                    false, //transposeA
                    *raw_type(P),
                    false, //transposeP
                    // *result.raw_type(),
                    *raw_result,
                    true  //call_FillComplete_on_result
                );

                result.raw_type() = raw_result;

            } catch(const std::exception &ex) {
                std::cerr << "RAP: " << ex.what() << std::endl;
                assert(false);
            }

            UTOPIA_TRACE_END(expr);
        }
    };

}

#endif //WITH_TRILINOS_TPETRAEXT

#endif //UTOPIA_TRILINOS_EVAL_RAP_HPP

