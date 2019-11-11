#ifndef UTOPIA_INLINE_EVAL_HPP
#define UTOPIA_INLINE_EVAL_HPP

#include "utopia_Traits.hpp"

namespace utopia {

    template<class Left, class Right, int Order = Traits<Left>::Order>
    class InlineEval {};

    template<class Left, class Right>
    class InlineEval<Left, Right, 1> {
    public:
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static void apply(Left &left, const Right &right)
        {
            const auto n = left.size();
            for(SizeType i = 0; i < n; ++i) {
                left(i) = right(i);
            }
        }
    };

    template<class Left, class Right>
    class InlineEval<Left, Right, 2> {
    public:
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static void apply(Left &left, const Right &right)
        {
            const auto rows = left.rows();
            const auto cols = left.cols();

            for(SizeType i = 0; i < rows; ++i) {
                for(SizeType j = 0; j < cols; ++j) {
                    left(i, j) = right(i, j);
                }
            }
        }
    };
}

#endif //UTOPIA_INLINE_EVAL_HPP
