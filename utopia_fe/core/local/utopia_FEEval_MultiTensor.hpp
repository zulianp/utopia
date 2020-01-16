#ifndef UTOPIA_FE_EVAL_MULTI_TENSOR_HPP
#define UTOPIA_FE_EVAL_MULTI_TENSOR_HPP

#include "utopia_Eval_Determinant.hpp"
#include "utopia_Eval_Inverse.hpp"
#include "utopia_LocalFEForwardDeclarations.hpp"

namespace utopia {

    template<class T, int Backend>
    class EvalDeterminant<MultiTensor<T, 2>, Backend> {
    public:
        using SizeType = typename Traits<MultiTensor<T, 2>>::SizeType;

        template<class Result>
        inline static bool apply(const MultiTensor<T, 2> &m, Result &result)
        {
            const SizeType n = m.size();

            if(n != result.size()) {
                result.resize(n);
            }

            for(SizeType i = 0; i < n; ++i) {
                result[i] = det(m[i]);
            }

            return true;
        }
    };

    template<class T, class Traits, int Backend>
    class Eval<Determinant<
                    Tensor<
                        MultiTensor<T, 2>,
                        2
                    >
                >, Traits, Backend> {
    public:
        using MT = MultiTensor<T, 2>;
        using T2 = utopia::Tensor<MT, 2>;
        using Expr = utopia::Determinant<T2>;
        using Scalar = typename Traits::Scalar;

        template<class Result>
        inline static void apply(const Expr &expr, Result &num)
        {
            UTOPIA_TRACE_BEGIN(expr);
            EvalDeterminant<MT>::apply(expr.expr().derived(), num);
            UTOPIA_TRACE_END(expr);
        }
    };

    template<class T, int Backend>
    class EvalInverse<MultiTensor<T, 2>, Backend> {
    public:
        using SizeType = typename Traits<MultiTensor<T, 2>>::SizeType;

        inline static bool apply(const MultiTensor<T, 2> &m, MultiTensor<T, 2> &result)
        {
            const SizeType n = m.size();

            if(n != result.size()) {
                result.resize(n);
            }

            for(SizeType i = 0; i < n; ++i) {
                result[i] = inv(m[i]);
            }

            return true;
        }
    };

}

#endif //UTOPIA_FE_EVAL_MULTI_TENSOR_HPP
