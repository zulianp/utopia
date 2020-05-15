#ifndef UTOPIA_FE_EVAL_REDUCE_HPP
#define UTOPIA_FE_EVAL_REDUCE_HPP

#include "utopia_AssemblyContext.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_FEBackend.hpp"
#include "utopia_FindSpace.hpp"

namespace utopia {

    template <class Traits, int Order>
    class InnerProduct {};

    template <class Traits>
    class InnerProduct<Traits, 0> {
    public:
        typedef typename Traits::Scalar ScalarType;
        typedef std::vector<ScalarType> Type;
        typedef typename Traits::Implementation Space;
        static const int Backend = Traits::Backend;

        template <class Left, class Right>
        inline static std::vector<ScalarType> apply(const Left &left,
                                                    const Right &right,
                                                    AssemblyContext<Traits::Backend> &ctx) {
            return FEBackend<Traits::Backend>::inner(FEEval<Left, Traits, Backend, QUAD_DATA_NO>::apply(left, ctx),
                                                     FEEval<Right, Traits, Backend, QUAD_DATA_NO>::apply(right, ctx),
                                                     ctx);
        }

        // inline static const Type & collect(const Type &out)
        // {
        // 	return out;
        // }

        // inline static Type collect(const std::vector<Type> &out)
        // {
        // 	Type reduced_out = out[0];

        // 	const std::size_t n = out.size();
        // 	for(std::size_t i = 1; i < n; ++i) {
        // 		reduced_out += out[i];
        // 	}

        // 	return reduced_out;
        // }
    };

    template <class Traits>
    class InnerProduct<Traits, 1> {
    public:
        typedef typename Traits::Vector Type;
        typedef typename Traits::Implementation Space;
        static const int Backend = Traits::Backend;

        template <class LeftExpr, class Test>
        inline static Type apply(const LeftExpr &left, const Test &test, AssemblyContext<Traits::Backend> &ctx) {
            auto of = offset(test, ctx);
            return FEBackend<Traits::Backend>::linear_form(
                of,
                FEEval<LeftExpr, Traits, Backend, QUAD_DATA_YES>::apply(left, ctx),
                FEEval<Test, Traits, Backend, QUAD_DATA_YES>::apply(test, ctx),
                ctx);
        }

        template <class Test>
        inline static Range offset(const Test &test, const AssemblyContext<Backend> &ctx) {
            FindTestSpace<Space, Test> f_test;
            f_test.apply(test);

            unsigned int n_shape_functions = 0;
            if (f_test.prod_space()) {
                n_shape_functions = FEBackend<Traits::Backend>::n_shape_functions(*f_test.prod_space(), ctx);
            } else {
                n_shape_functions = FEBackend<Traits::Backend>::n_shape_functions(*f_test.space(), ctx);
            }

            const unsigned int offset = FEBackend<Traits::Backend>::offset(*f_test.space(), ctx);
            return Range(offset, offset + n_shape_functions);
        }
    };

    template <class Traits>
    class InnerProduct<Traits, 2> {
    public:
        typedef typename Traits::Matrix Type;
        typedef typename Traits::Implementation Space;
        static const int Backend = Traits::Backend;

        template <class Trial, class Test>
        inline static Type apply(const Trial &trial, const Test &test, AssemblyContext<Traits::Backend> &ctx) {
            auto of = offsets(trial, test, ctx);
            return FEBackend<Traits::Backend>::bilinear_form(
                of.first,
                of.second,
                FEEval<Trial, Traits, Backend, QUAD_DATA_YES>::apply(trial, ctx),
                FEEval<Test, Traits, Backend, QUAD_DATA_YES>::apply(test, ctx),
                ctx);
        }

        template <class Trial, class Test>
        inline static std::pair<Range, Range> offsets(const Trial &trial,
                                                      const Test &test,
                                                      const AssemblyContext<Backend> &ctx) {
            FindTestSpace<Space, Test> f_test;
            f_test.apply(test);

            unsigned int n_test_functions = 0;
            if (f_test.prod_space()) {
                n_test_functions = FEBackend<Traits::Backend>::n_shape_functions(*f_test.prod_space(), ctx);
            } else {
                n_test_functions = FEBackend<Traits::Backend>::n_shape_functions(*f_test.space(), ctx);
            }

            const unsigned int test_offset = FEBackend<Traits::Backend>::offset(*f_test.space(), ctx);

            FindTrialSpace<Space, Trial> f_trial;
            f_trial.apply(trial);

            unsigned int n_trial_functions = 0;
            if (f_trial.prod_space()) {
                n_trial_functions = FEBackend<Traits::Backend>::n_shape_functions(*f_trial.prod_space(), ctx);
            } else {
                n_trial_functions = FEBackend<Traits::Backend>::n_shape_functions(*f_trial.space(), ctx);
            }

            const unsigned int trial_offset = FEBackend<Traits::Backend>::offset(*f_trial.space(), ctx);

            return {Range(trial_offset, trial_offset + n_trial_functions),
                    Range(test_offset, test_offset + n_test_functions)};
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class FEEval<Reduce<Binary<Left, Right, EMultiplies>, Plus>, Traits, Backend, QUAD_DATA_YES> {
    public:
        typedef utopia::Reduce<Binary<Left, Right, EMultiplies>, Plus> Expr;

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Traits::Backend>::inner(
                FEEval<Left, Traits, Backend, QUAD_DATA_YES>::apply(expr.expr().left(), ctx),
                FEEval<Right, Traits, Backend, QUAD_DATA_YES>::apply(expr.expr().right(), ctx),
                ctx)) {
            return FEBackend<Traits::Backend>::inner(
                FEEval<Left, Traits, Backend, QUAD_DATA_YES>::apply(expr.expr().left(), ctx),
                FEEval<Right, Traits, Backend, QUAD_DATA_YES>::apply(expr.expr().right(), ctx),
                ctx);
        }
    };

    template <class Left, class Right, class Traits, int Backend>
    class FEEval<Reduce<Binary<Left, Right, EMultiplies>, Plus>, Traits, Backend, QUAD_DATA_NO> {
    public:
        typedef utopia::Reduce<Binary<Left, Right, EMultiplies>, Plus> Expr;
        typedef utopia::IsBilinearFormSplit<Left, Right> IsSplit;

        typedef utopia::InnerProduct<Traits, IsSplit::order> InnerProductT;
        typedef typename InnerProductT::Type Result;

        static_assert(IsSplit::order == 0 || IsSplit::has_test, "If it is a form it must have a test function");

        inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx) -> Result {
            // For handling multilinear forms:
            //- is it a form or just a reduce? if reduce eval and return
            //- identify if it is a simple form (use ctx to do that) and handle simple case first
            //- if it is part of a composite equation then find out which "sub-block" (not necessarily contiguous)
            //  is this form about. create index sets and write entries in block (implement in backend only)

            if (IsSplit::right_has_test) {
                return InnerProductT::apply(expr.expr().left(), expr.expr().right(), ctx);
            } else {
                return InnerProductT::apply(expr.expr().right(), expr.expr().left(), ctx);
            }
        }
    };

    template <class Expr, class AssemblyContext>
    class FunctionalTraits<Reduce<Expr, Plus>, AssemblyContext> {
    public:
        inline static int type(const Reduce<Expr, Plus> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Expr, AssemblyContext>::type(expr.expr(), ctx);
        }

        inline static int order(const Reduce<Expr, Plus> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Expr, AssemblyContext>::order(expr.expr(), ctx);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_REDUCE_HPP
