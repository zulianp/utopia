#ifndef UTOPIA_VIEW_ASSIGN_HPP
#define UTOPIA_VIEW_ASSIGN_HPP

#include "utopia_Base.hpp"
#include "utopia_DeviceExpression.hpp"
#include "utopia_InlineEval.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class Left, class Right>
    class DeviceAssign {
    public:
        UTOPIA_INLINE_FUNCTION DeviceAssign(Left &left, const Right &right) : left_(left), right_(right) {}

        UTOPIA_INLINE_FUNCTION void apply() { apply(left_, right_); }

        UTOPIA_INLINE_FUNCTION static void apply(Left &left, const Right &right) {
            InlineEval<Left, Right>::apply(left, right);
        }

    private:
        Left &left_;
        UTOPIA_STORE_CONST(Right) right_;
    };

    template <class Left, class Right>
    class DeviceAssign<Left, DeviceInverse<Right>> {
    public:
        UTOPIA_INLINE_FUNCTION static void apply(Left &left, const DeviceInverse<Right> &right) { right.apply(left); }
    };

}  // namespace utopia

#endif