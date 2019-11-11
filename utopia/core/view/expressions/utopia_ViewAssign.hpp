#ifndef UTOPIA_VIEW_ASSIGN_HPP
#define UTOPIA_VIEW_ASSIGN_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_InlineEval.hpp"
#include "utopia_ViewExpression.hpp"

namespace utopia {

    template<class Left, class Right>
    class ViewAssign {
    public:
        UTOPIA_INLINE_FUNCTION ViewAssign(Left &left, const Right &right) : left_(left), right_(right) {}

        UTOPIA_INLINE_FUNCTION void apply()
        {
            InlineEval<Left, Right>::apply(left_, right_);
        }

    private:
        Left &left_;
        UTOPIA_STORE_CONST(Right) right_;

    };
}

#endif