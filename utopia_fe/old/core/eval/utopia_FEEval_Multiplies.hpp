#ifndef UTOPIA_FE_EVAL_MULTIPLIES_HPP
#define UTOPIA_FE_EVAL_MULTIPLIES_HPP

#include "utopia_AssemblyContext.hpp"
#include "utopia_Eval_Empty.hpp"
#include "utopia_FEBackend.hpp"
#include "utopia_FEEval_Binary.hpp"
#include "utopia_FEEval_Reduce.hpp"

namespace utopia {
    template <class Left, class Right, class Traits, int Backend, int IsHackType>
    class MultipliesDelegate {
    public:
        inline static auto apply(const Left &left, const Right &right, AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend, IsHackType>::apply(left, ctx),
                                                         FEEval<Right, Traits, Backend, IsHackType>::apply(right, ctx),
                                                         Multiplies(),
                                                         ctx)) {
            return FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend, IsHackType>::apply(left, ctx),
                                                    FEEval<Right, Traits, Backend, IsHackType>::apply(right, ctx),
                                                    Multiplies(),
                                                    ctx);
        }
    };

    template <typename T, class Left, class Right, class Traits, int Backend, int IsQuadData>
    class MultipliesDelegate<Number<T>,
                             Binary<Transposed<Gradient<Left> >, Gradient<Right>, Plus>,
                             Traits,
                             Backend,
                             IsQuadData> {
    public:
        typedef utopia::Binary<Transposed<Gradient<Left> >, Gradient<Right>, Plus> SymmetrizedGrad;

        inline static auto apply(const Number<T> &scale_factor,
                                 const SymmetrizedGrad &symm_grad,
                                 AssemblyContext<Backend> &ctx)
            -> decltype(FEBackend<Backend>::grad_t_plus_grad(scale_factor,
                                                             symm_grad.left().expr().expr(),
                                                             symm_grad.right().expr(),
                                                             ctx)) {
            return FEBackend<Backend>::grad_t_plus_grad(
                scale_factor, symm_grad.left().expr().expr(), symm_grad.right().expr(), ctx);
        }
    };

    template <class Left, class Right, class Traits, int Backend, int IsQuadData>
    class MultipliesDelegate<Number<Left>, Integral<Right>, Traits, Backend, IsQuadData> {
    public:
        inline static auto apply(const Number<Left> &left, const Integral<Right> &right, AssemblyContext<Backend> &ctx)
            -> decltype(FEEval<Right, Traits, Backend, IsQuadData>::apply(right.expr(), ctx)) {
            return static_cast<Left>(left) * FEEval<Right, Traits, Backend, IsQuadData>::apply(right.expr(), ctx);
        }
    };

    // template<class Left, class Right, class Op, class Traits, int Backend, int IsQuadData>
    // class MultipliesDelegate<
    // 			Number<Left>,
    // 			Reduce<Right, Op>,
    // 			Traits,
    // 			Backend,
    // 			IsQuadData
    // 		> {
    // public:
    //     inline static auto apply(
    //     	const Number<Left> &left,
    //     	const Reduce<Right, Op> &right,
    //     	AssemblyContext<Backend> &ctx) -> decltype(
    //     	    		FEEval<Reduce<Right, Op>, Traits, Backend, IsQuadData>::apply(right, ctx)
    //     	    	)
    //     {
    //     	return static_cast<Left>(left) * FEEval<Reduce<Right, Op>, Traits, Backend, IsQuadData>::apply(right,
    //     ctx);
    //     }
    // };

    template <class Left, class Right, class Traits, int Backend, int IsQuadData>
    class MultipliesDelegate<Integral<Left>, Number<Right>, Traits, Backend, IsQuadData> {
    public:
        inline static auto apply(const Integral<Left> &left, const Number<Right> &right, AssemblyContext<Backend> &ctx)
            -> decltype(FEEval<Left, Traits, Backend, IsQuadData>::apply(left, ctx)) {
            return FEEval<Left, Traits, Backend, IsQuadData>::apply(left, ctx) * static_cast<Right>(right);
        }
    };

    template <class T, class C, class F, class Traits, int Backend, int IsQuadData>
    class MultipliesDelegate<Interpolate<C, F>, ConstantCoefficient<T, 0>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::Interpolate<C, F> Left;
        typedef utopia::ConstantCoefficient<T, 0> Right;
        typedef utopia::Binary<Left, Right, Multiplies> Expr;

        static auto apply(const Left &left, const Right &right, AssemblyContext<Backend> &ctx)
            -> decltype(FEEval<Left, Traits, Backend, IsQuadData>::apply(left, ctx)) {
            return FEBackend<Backend>::apply_binary(right, left, Multiplies(), ctx);
        }
    };

    // see Reduce
    template <class Left, class Right, class Traits, int Backend>
    class MultipliesDelegate<Left, Right, Traits, Backend, -1> {
    public:
        typedef utopia::IsBilinearFormSplit<Left, Right> IsSplit;

        typedef utopia::InnerProduct<Traits, IsSplit::order> InnerProductT;
        typedef typename InnerProductT::Type Result;

        inline static auto apply(const Left &left, const Right &right, AssemblyContext<Backend> &ctx) -> Result {
            if (IsSplit::right_has_test) {
                return InnerProductT::apply(left, right, ctx);
            } else {
                return InnerProductT::apply(right, left, ctx);
            }
        }
    };

    template <class Left, class Right, class Traits, int Backend, int IsQuadData>
    class FEEval<Binary<Left, Right, Multiplies>, Traits, Backend, IsQuadData> {
    public:
        typedef utopia::IsBilinearFormSplit<Left, Right> IsSplit;
        static const int IsHackType = (IsSplit::value) ? -1 : IsQuadData;

        typedef utopia::MultipliesDelegate<Left, Right, Traits, Backend, IsHackType> DelegateT;

        inline static auto apply(const Binary<Left, Right, Multiplies> &expr, AssemblyContext<Backend> &ctx)
            -> decltype(DelegateT::apply(expr.left(), expr.right(), ctx)) {
            return DelegateT::apply(expr.left(), expr.right(), ctx);
        }
    };

    template <class Left, class Right, class AssemblyContext>
    class FunctionalTraits<Binary<Left, Right, Multiplies>, AssemblyContext> {
    public:
        inline static int type(const Binary<Left, Right, Multiplies> &expr, const AssemblyContext &ctx) {
            return std::max(FunctionalTraits<Left, AssemblyContext>::type(expr.left(), ctx),
                            FunctionalTraits<Right, AssemblyContext>::type(expr.right(), ctx));
        }

        inline static int order(const Binary<Left, Right, Multiplies> &expr, const AssemblyContext &ctx) {
            return FunctionalTraits<Left, AssemblyContext>::order(expr.left(), ctx) +
                   FunctionalTraits<Right, AssemblyContext>::order(expr.right(), ctx);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_MULTIPLIES_HPP
