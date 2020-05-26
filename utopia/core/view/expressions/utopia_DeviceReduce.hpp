#ifndef UTOPIA_DEVICE_REDUCE_HPP
#define UTOPIA_DEVICE_REDUCE_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_Base.hpp"
#include "utopia_DeviceOp.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Expr, class Op, int Order = Traits<Expr>::Order>
    class DeviceReduce;

    template <class Expr, class Op>
    class DeviceReduce<Expr, Op, 1> {
    public:
        using Scalar = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &expr, const Scalar &initial_value) {
            Scalar ret = initial_value;

            const SizeType n = expr.size();

            for (SizeType i = 0; i < n; ++i) {
                ret = DeviceOp<Scalar, Op>::apply(ret, expr(i));
            }

            return ret;
        }

        UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &expr) {
            const SizeType n = expr.size();
            Scalar ret = expr(0);

            for (SizeType i = 1; i < n; ++i) {
                ret = DeviceOp<Scalar, Op>::apply(ret, expr(i));
            }

            return ret;
        }
    };

    template <class Expr, class Op>
    class DeviceReduce<Expr, Op, 2> {
    public:
        using Scalar = typename Traits<Expr>::Scalar;
        using SizeType = typename Traits<Expr>::SizeType;

        UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &expr) {
            const SizeType rows = expr.rows();
            const SizeType cols = expr.cols();

            Scalar ret = expr(0, 0);

            for (SizeType j = 1; j < cols; ++j) {
                ret = DeviceOp<Scalar, Op>::apply(ret, expr(0, j));
            }

            for (SizeType i = 1; i < rows; ++i) {
                for (SizeType j = 0; j < cols; ++j) {
                    ret = DeviceOp<Scalar, Op>::apply(ret, expr(i, j));
                }
            }

            return ret;
        }

        UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &expr, const Scalar &initial_value) {
            Scalar ret = initial_value;

            const SizeType rows = expr.rows();
            const SizeType cols = expr.cols();

            for (SizeType i = 0; i < rows; ++i) {
                for (SizeType j = 0; j < cols; ++j) {
                    ret = DeviceOp<Scalar, Op>::apply(ret, expr(i, j));
                }
            }

            return ret;
        }
    };

    // template<class Expr, class Op>
    // class DeviceReduce<Expr, Op, 4> {
    // public:
    //     using Scalar   = typename Traits<Expr>::Scalar;
    //     using SizeType = typename Traits<Expr>::SizeType;

    //     UTOPIA_INLINE_FUNCTION static Scalar apply(const Expr &expr)
    //     {
    //         const SizeType N0 = extent(expr, 0);
    //         const SizeType N1 = extent(expr, 1);
    //         const SizeType N2 = extent(expr, 2);
    //         const SizeType N3 = extent(expr, 3);

    //         Scalar ret = expr(0, 0, 0, 0);

    //         for(SizeType j = 0; j < N1; ++j) {
    //             for(SizeType k = 0; k < N2; ++k) {
    //                 for(SizeType l = 1; l < N3; ++l) {
    //                     ret = DeviceOp<Scalar, Op>::apply(ret, expr(0, j, k, l));
    //                 }
    //             }
    //         }

    //         for(SizeType i = 1; i < rows; ++i) {
    //             for(SizeType j = 0; j < cols; ++j) {
    //                 ret = DeviceOp<Scalar, Op>::apply(ret, expr(i, j));
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
    //                 ret = DeviceOp<Scalar, Op>::apply(ret, expr(i, j));
    //             }
    //         }

    //         return ret;
    //     }
    // };

}  // namespace utopia

#endif