#ifndef UTOPIA_DEVICE_APPROX_EQUAL_HPP
#define UTOPIA_DEVICE_APPROX_EQUAL_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    template <class Left, class Right, int Order = Traits<Left>::Order>
    class DeviceApproxEqual {};

    template <class Left, class Right>
    class DeviceApproxEqual<Left, Right, 1> {
    public:
        using Scalar = typename Traits<Left>::Scalar;
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static bool apply(const Left &left, const Right &right, const Scalar &tol) {
            const SizeType n = left.size();
            bool ret = n != 0;
            for (SizeType i = 0; i < n; ++i) {
                if (!device::approxeq(left(i), right(i), tol)) {
                    ret = false;
                    break;
                }
            }

            return ret;
        }
    };

    template <class Left, class Right>
    class DeviceApproxEqual<Left, Right, 2> {
    public:
        using Scalar = typename Traits<Left>::Scalar;
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static bool apply(const Left &left, const Right &right, const Scalar &tol) {
            const SizeType rows = left.rows();
            const SizeType cols = left.cols();

            bool ret = rows != 0;

            for (SizeType i = 0; i < rows; ++i) {
                for (SizeType j = 0; j < cols; ++j) {
                    if (!device::approxeq(left(i, j), right(i, j), tol)) {
                        ret = false;
                        break;
                    }
                }
            }

            return ret;
        }
    };

    template <class Left, class Right>
    class DeviceApproxEqual<Left, Right, 4> {
    public:
        using Scalar = typename Traits<Left>::Scalar;
        using SizeType = typename Traits<Left>::SizeType;

        UTOPIA_INLINE_FUNCTION static bool apply(const Left &left, const Right &right, const Scalar &tol) {
            const SizeType N0 = device::extent(left, 0);
            const SizeType N1 = device::extent(left, 1);
            const SizeType N2 = device::extent(left, 2);
            const SizeType N3 = device::extent(left, 3);

            if (N0 != device::extent(right, 0) || N1 != device::extent(right, 1) || N2 != device::extent(right, 2) ||
                N3 != device::extent(right, 3)) {
                return false;
            }

            bool ret = true;
            for (SizeType i = 0; i < N0; ++i) {
                for (SizeType j = 0; j < N1; ++j) {
                    for (SizeType k = 0; k < N2; ++k) {
                        for (SizeType l = 0; l < N3; ++l) {
                            if (!device::approxeq(left(i, j, k, l), right(i, j, k, l), tol)) {
                                ret = false;
                                break;
                            }
                        }
                    }
                }
            }

            return ret;
        }
    };

    template <class Left, class Right>
    class DeviceApproxEqual<Left, Right, 0> {
    public:
        using Scalar = typename Traits<Left>::Scalar;

        UTOPIA_INLINE_FUNCTION bool apply(const Left &left, const Right &right, const Scalar &tol) {
            return device::approxeq(left, right, tol);
        }
    };

    template <class Left, class Right>
    class Traits<DeviceApproxEqual<Left, Right> > : public Traits<typename MostDescriptive<Left, Right>::Type> {};
}  // namespace utopia

#endif
