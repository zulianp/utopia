#ifndef UTOPIA_DEVICE_DIAG_HPP
#define UTOPIA_DEVICE_DIAG_HPP

#include "utopia_Base.hpp"
#include "utopia_DeviceExpression.hpp"
#include "utopia_DeviceOp.hpp"
#include "utopia_InlineEval.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <class InnerExpr>
    class DeviceDiag : public DeviceExpression<DeviceDiag<InnerExpr>> {
    public:
        using SizeType = typename Traits<InnerExpr>::SizeType;
        using Scalar = typename Traits<InnerExpr>::Scalar;

        UTOPIA_INLINE_FUNCTION DeviceDiag(const InnerExpr &expr) : expr_(expr) {}

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const {
            return (i == j) ? expr_(i) : 0.0;
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i) const { return expr_(i, i); }

        UTOPIA_INLINE_FUNCTION SizeType rows() const { return expr_.size(); }

        UTOPIA_INLINE_FUNCTION SizeType cols() const { return expr_.size(); }

        UTOPIA_INLINE_FUNCTION SizeType size() const { return device::min(expr_.rows(), expr_.cols()); }

    private:
        UTOPIA_STORE_CONST(InnerExpr) expr_;
    };

    template <class InnerExpr>
    class Traits<DeviceDiag<InnerExpr>> : public Traits<InnerExpr> {
    public:
        static const int Order = Traits<InnerExpr>::Order == 1 ? 2 : 1;
    };
}  // namespace utopia

#endif  // UTOPIA_DEVICEDIAG_HPP