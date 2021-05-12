#ifndef UTOPIA_DEVICE_MULTIPLY_HPP
#define UTOPIA_DEVICE_MULTIPLY_HPP

#include "utopia_Base.hpp"

#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ViewForwardDeclarations.hpp"

#include "utopia_Algorithms.hpp"

namespace utopia {

    template <class Left, class Right>
    class DeviceMultiply : public DeviceExpression<DeviceMultiply<Left, Right>> {
    public:
        using SizeType = typename Traits<Left>::SizeType;
        using Scalar = typename std::remove_const<typename Traits<Left>::Scalar>::type;

        UTOPIA_INLINE_FUNCTION DeviceMultiply(const Left &left, const Right &right) : left_(left), right_(right) {
            // UTOPIA_DEVICE_ASSERT(left.cols() == right.size());
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i, const SizeType &j) const {
            const SizeType cols = left_.cols();

            Scalar ret = 0.0;
            for (SizeType k = 0; k < cols; ++k) {
                ret += left_(i, k) * right_(k, j);
            }

            return ret;
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i,
                                                 const SizeType &j,
                                                 const SizeType &k,
                                                 const SizeType &l) const {
            // http://www.mate.tue.nl/~peters/4K400/VectTensColMat.pdf
            assert(device::extent(left_, 0) == device::extent(left_, 1));
            assert(device::extent(left_, 0) == device::extent(left_, 2));
            assert(device::extent(left_, 0) == device::extent(left_, 3));
            assert(device::extent(left_, 0) == device::extent(right_, 0));
            assert(device::extent(left_, 1) == device::extent(right_, 1));
            assert(device::extent(left_, 2) == device::extent(right_, 2));
            assert(device::extent(left_, 3) == device::extent(right_, 3));

            const SizeType N = device::extent(left_, 0);

            Scalar sum = 0.0;
            for (SizeType m = 0; m < N; ++m) {
                for (SizeType n = 0; n < N; ++n) {
                    sum += left_(i, j, m, n) * right_(n, m, k, l);
                }
            }

            return sum;
        }

        UTOPIA_INLINE_FUNCTION Scalar operator()(const SizeType &i) const {
            const SizeType cols = left_.cols();

            Scalar ret = 0.0;
            for (SizeType j = 0; j < cols; ++j) {
                ret += left_(i, j) * right_(j);
            }

            return ret;
        }

        inline std::string get_class() const override {
            return "DeviceMultiply<" + left_.get_class() + ", " + right_.get_class() + ">";
        }

        UTOPIA_INLINE_FUNCTION SizeType rows() const { return left_.rows(); }

        UTOPIA_INLINE_FUNCTION SizeType cols() const { return right_.cols(); }

    private:
        UTOPIA_STORE_CONST(Left) left_;
        UTOPIA_STORE_CONST(Right) right_;
    };

    template <class Left, class Right>
    class Traits<DeviceMultiply<Left, Right>> : public Traits<typename ChooseType<Left, Right, Right>::Type> {};
}  // namespace utopia

#endif
