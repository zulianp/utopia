#ifndef UTOPIA_INLINE_EVAL_HPP
#define UTOPIA_INLINE_EVAL_HPP

#include "utopia_Traits.hpp"

namespace utopia {

    template <class Left, class Right, int Order = Traits<Left>::Order>
    class InlineEval {};

    template <class Left, class Right>
    class InlineEval<Left, Right, 1> {
    public:
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static void apply(Left &left, const Right &right) {
            const SizeType n = left.size();
            for (SizeType i = 0; i < n; ++i) {
                left(i) = right(i);
            }
        }
    };

    template <class Left, class Right>
    class InlineEval<Left, Right, 2> {
    public:
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static void apply(Left &left, const Right &right) {
            const SizeType rows = left.rows();
            const SizeType cols = left.cols();

            for (SizeType i = 0; i < rows; ++i) {
                for (SizeType j = 0; j < cols; ++j) {
                    left(i, j) = right(i, j);
                }
            }
        }
    };

    template <class Left, class Right>
    class InlineEval<Left, Right, 4> {
    public:
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static void apply(Left &left, const Right &right) {
            const SizeType N0 = extent(left, 0);
            const SizeType N1 = extent(left, 1);
            const SizeType N2 = extent(left, 2);
            const SizeType N3 = extent(left, 3);

            for (SizeType i = 0; i < N0; ++i) {
                for (SizeType j = 0; j < N1; ++j) {
                    for (SizeType k = 0; k < N2; ++k) {
                        for (SizeType l = 0; l < N3; ++l) {
                            left(i, j, k, l) = right(i, j, k, l);
                        }
                    }
                }
            }
        }
    };
}  // namespace utopia

#endif  // UTOPIA_INLINE_EVAL_HPP
