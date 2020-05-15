#ifndef UTOPIA_DEVICE_REDUCE_FIND_HPP
#define UTOPIA_DEVICE_REDUCE_FIND_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Expr, int Order = Traits<Expr>::Order>
    class DeviceFind;

    template <class Expr>
    class DeviceFind<Expr, 1> {
    public:
        using Scalar = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        template <class Predicate>
        UTOPIA_INLINE_FUNCTION static void apply(const Expr &expr,
                                                 Predicate pred,
                                                 // OUTPUT
                                                 SizeType &idx,
                                                 Scalar &value) {
            const SizeType n = expr.size();

            value = expr(0);
            idx = 0;

            for (SizeType i = 1; i < n; ++i) {
                const Scalar val_i = expr(i);
                if (pred(val_i, value)) {
                    value = val_i;
                    idx = i;
                }
            }
        }
    };

    template <class Expr>
    typename DeviceFind<Expr, 1>::SizeType imax(const DeviceExpression<Expr> &expr) {
        using Scalar = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        Scalar val = 0.0;
        SizeType idx = 0;

        DeviceFind<Expr, 1>::apply(
            expr.derived(), UTOPIA_LAMBDA(const Scalar &l, const Scalar &r)->bool { return l > r; }, idx, val);

        return idx;
    }

    // template<class Expr, class Predicate>
    // class DeviceFind<Expr, Predicate, 2> {
    // public:
    //     using Scalar   = typename Traits<Expr>::Scalar;
    //     using SizeType = typename Traits<Expr>::SizeType;

    //     UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &expr)
    //     {
    //         const SizeType rows = expr.rows();
    //         const SizeType cols = expr.cols();

    //         Scalar ret = expr(0, 0);
    //         SizeType idx =

    //         for(SizeType j = 1; j < cols; ++j) {
    //             ret = DevicePredicate<Scalar, Predicate>::apply(ret, expr(0, j));
    //         }

    //         for(SizeType i = 1; i < rows; ++i) {
    //             for(SizeType j = 0; j < cols; ++j) {
    //                 ret = DevicePredicate<Scalar, Predicate>::apply(ret, expr(i, j));
    //             }
    //         }

    //         return ret;
    //     }

    //     UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &expr, const Scalar &initial_value)
    //     {
    //         Scalar ret = initial_value;

    //         const SizeType rows = expr.rows();
    //         const SizeType cols = expr.cols();

    //         for(SizeType i = 0; i < rows; ++i) {
    //             for(SizeType j = 0; j < cols; ++j) {
    //                 ret = DevicePredicate<Scalar, Predicate>::apply(ret, expr(i, j));
    //             }
    //         }

    //         return ret;
    //     }
    // };

}  // namespace utopia
#endif  // UTOPIA_DEVICE_REDUCE_FIND_HPP
