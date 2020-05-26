#ifndef UTOPIA_AUTO_DIFF_HPP
#define UTOPIA_AUTO_DIFF_HPP

#include "utopia_AutoDiffExpr.hpp"
#include "utopia_Differentiable.hpp"
#include "utopia_Simplify.hpp"

#include "utopia_AutoDiffExpr_Binary.hpp"
#include "utopia_AutoDiffExpr_Multiply.hpp"
#include "utopia_AutoDiffExpr_Reduce.hpp"
#include "utopia_AutoDiffExpr_Trace.hpp"
#include "utopia_AutoDiffExpr_Transposed.hpp"

#include "utopia_Simplify_Binary.hpp"
#include "utopia_Simplify_Multiply.hpp"

namespace utopia {
    template <class Derived>
    auto simplify(const Expression<Derived> &expr) -> decltype(Simplify<Derived>::make(expr.derived())) {
        return Simplify<Derived>::make(expr.derived());
    }
}  // namespace utopia

#endif  // UTOPIA_AUTO_DIFF_HPP
