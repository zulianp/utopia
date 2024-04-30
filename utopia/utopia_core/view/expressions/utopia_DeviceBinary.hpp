#ifndef UTOPIA_VIEW_BINARY_HPP
#define UTOPIA_VIEW_BINARY_HPP

#include "utopia_Base.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ViewForwardDeclarations.hpp"

#include "utopia_DeviceExpression.hpp"
#include "utopia_DeviceOp.hpp"
#include "utopia_InlineEval.hpp"

namespace utopia {

    template <class T>
    class HasSize {
    public:
        static const int value = 1;
    };

    template <class Left, class Right, int LeftHasSize = HasSize<Left>::value, int RightHasSize = HasSize<Right>::value>
    class GetSize {
    public:
        using SizeType = typename Traits<Right>::SizeType;

        UTOPIA_INLINE_FUNCTION static constexpr SizeType rows(const Left &, const Right &r) { return r.rows(); }

        UTOPIA_INLINE_FUNCTION static constexpr SizeType cols(const Left &, const Right &r) { return r.cols(); }

        UTOPIA_INLINE_FUNCTION static constexpr SizeType extent(const Left &, const Right &r, const SizeType i) {
            return device::extent(r, i);
        }
    };

    template <class Left, class Right>
    class GetSize<Left, Right, 1, 0> {
    public:
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static constexpr SizeType rows(const Left &l, const Right &) { return l.rows(); }

        UTOPIA_INLINE_FUNCTION static constexpr SizeType cols(const Left &l, const Right &) { return l.cols(); }

        UTOPIA_INLINE_FUNCTION static constexpr SizeType extent(const Left &l, const Right &, const SizeType i) {
            return device::extent(l, i);
        }
    };

    template <class Left, class Right, class Op>
    class DeviceBinary : public DeviceExpression<DeviceBinary<Left, Right, Op>> {
    public:
        using SizeType = typename Traits<Left>::SizeType;
        using Scalar = decltype(typename Traits<Left>::Scalar(0) + typename Traits<Right>::Scalar(0));

        UTOPIA_INLINE_FUNCTION constexpr DeviceBinary(const Left &left, const Right &right)
            : left_(left), right_(right) {
            // UTOPIA_DEVICE_ASSERT(left.size() == right.size());
        }

        UTOPIA_INLINE_FUNCTION constexpr Scalar operator()(const SizeType &i,
                                                           const SizeType &j,
                                                           const SizeType &k,
                                                           const SizeType &l) const

        {
            return DeviceOp<Scalar, Op>::apply(left_(i, j, k, l), right_(i, j, k, l));
        }

        UTOPIA_INLINE_FUNCTION constexpr Scalar operator()(const SizeType &i, const SizeType &j) const {
            return DeviceOp<Scalar, Op>::apply(left_(i, j), right_(i, j));
        }

        UTOPIA_INLINE_FUNCTION constexpr Scalar operator()(const SizeType &i) const {
            return DeviceOp<Scalar, Op>::apply(left_(i), right_(i));
        }

        inline std::string get_class() const override {
            return std::string("DeviceBinary<") + left_.get_class() + ", " + right_.get_class() + ", " +
                   GetClass<Op>() + ">";
        }

        inline constexpr SizeType size() const { return left_.size(); }

        UTOPIA_INLINE_FUNCTION constexpr const Left &left() const { return left_; }

        UTOPIA_INLINE_FUNCTION constexpr const Right &right() const { return right_; }

        UTOPIA_INLINE_FUNCTION constexpr SizeType rows() const { return GetSize<Left, Right>::rows(left_, right_); }

        UTOPIA_INLINE_FUNCTION constexpr SizeType cols() const { return GetSize<Left, Right>::cols(left_, right_); }

        UTOPIA_INLINE_FUNCTION constexpr SizeType extent(const int i) const {
            return GetSize<Left, Right>::extent(left_, right_, i);
        }

    private:
        UTOPIA_STORE_CONST(Left) left_;
        UTOPIA_STORE_CONST(Right) right_;
    };

    template <class Left, class Right, class Op>
    class DeviceBinary<DeviceNumber<Left>, Right, Op>
        : public DeviceExpression<DeviceBinary<DeviceNumber<Left>, Right, Op>> {
    public:
        using SizeType = typename Traits<Right>::SizeType;
        using Scalar = decltype(typename Traits<Right>::Scalar(0) + Left(0));

        UTOPIA_INLINE_FUNCTION DeviceBinary(const DeviceNumber<Left> &left, const Right &right)
            : left_(left), right_(right) {}

        UTOPIA_INLINE_FUNCTION Scalar
        operator()(const SizeType &i, const SizeType &j, const SizeType &k, const SizeType &l) const

        {
            return DeviceOp<Scalar, Op>::apply(left_, right_(i, j, k, l));
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const {
            return DeviceOp<Scalar, Op>::apply(left_, right_(i, j));
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i) const {
            return DeviceOp<Scalar, Op>::apply(left_, right_(i));
        }

        UTOPIA_INLINE_FUNCTION const Left &left() const { return left_; }

        UTOPIA_INLINE_FUNCTION const Right &right() const { return right_; }

        UTOPIA_INLINE_FUNCTION SizeType rows() const { return right_.rows(); }

        UTOPIA_INLINE_FUNCTION SizeType cols() const { return right_.cols(); }

        UTOPIA_INLINE_FUNCTION constexpr SizeType extent(const int i) const { return device::extent(right_, i); }

    private:
        const Left left_;
        UTOPIA_STORE_CONST(Right) right_;
    };

    template <class Left, class Right, class Op>
    class HasSize<DeviceBinary<Left, Right, Op>> {
    public:
        static const int value = HasSize<Left>::value || HasSize<Right>::value;
    };

    template <class Left, class Right, class Operation>
    class Traits<DeviceBinary<Left, Right, Operation>> : public Traits<typename MostDescriptive<Left, Right>::Type> {};

    template <class Left, class Right, class Operation>
    UTOPIA_INLINE_FUNCTION typename Traits<DeviceBinary<Left, Right, Operation>>::SizeType rows(
        const DeviceBinary<Left, Right, Operation> &expr) {
        return expr.rows();
    }

    template <class Left, class Right, class Operation>
    UTOPIA_INLINE_FUNCTION typename Traits<DeviceBinary<Left, Right, Operation>>::SizeType cols(
        const DeviceBinary<Left, Right, Operation> &expr) {
        return expr.cols();
    }
}  // namespace utopia

#endif
