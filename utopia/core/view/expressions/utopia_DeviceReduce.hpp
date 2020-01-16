#ifndef UTOPIA_DEVICE_REDUCE_HPP
#define UTOPIA_DEVICE_REDUCE_HPP

#include "utopia_Base.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template<class Expr, class Op, int Order = Traits<Expr>::Order>
    class DeviceReduce;

    template<class Expr, class Op>
    class DeviceReduce<Expr, Op, 1> {
    public:
        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &expr, const Scalar &initial_value)
        {
            Scalar ret = initial_value;

            const SizeType n = expr.size();

            for(SizeType i = 0; i < n; ++i) {
                ret += DeviceOp<Scalar, Op>::apply(ret, expr(i));
            }

            return ret;
        }
    };

    template<class Expr, class Op>
    class DeviceReduce<Expr, Op, 2> {
    public:
        using Scalar   = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &expr, const Scalar &initial_value)
        {
            Scalar ret = initial_value;

            const SizeType rows = expr.rows();
            const SizeType cols = expr.cols();

            for(SizeType i = 0; i < rows; ++i) {
                for(SizeType j = 0; j < cols; ++j) {
                    ret += DeviceOp<Scalar, Op>::apply(ret, expr(i, j));
                }
            }

            return ret;
        }
    };

}

#endif