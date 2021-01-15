#ifndef UTOPIA_UNIFORM_TRI3_HPP
#define UTOPIA_UNIFORM_TRI3_HPP

#include "utopia_CppMacros.hpp"
#include "utopia_DeviceNumber.hpp"
#include "utopia_Edge2.hpp"
#include "utopia_Elem.hpp"
#include "utopia_MemType.hpp"
#include "utopia_Views.hpp"

namespace utopia {

    class RefTri3 {
    public:
        template <typename Point>
        UTOPIA_INLINE_FUNCTION static UTOPIA_CONSTEXPR auto fun(const int i, const Point &p) ->
            typename Traits<Point>::Scalar {
            const auto x = p[0];
            const auto y = p[1];

            switch (i) {
                case 0: {
                    return 1 - x - y;
                }
                case 1: {
                    return x;
                }
                case 2: {
                    return y;
                }
                default: {
                    UTOPIA_DEVICE_ASSERT_CXX14(false);
                    return 0.0;
                }
            }

            return 0.0;
        }

        template <typename Point>
        UTOPIA_INLINE_FUNCTION static UTOPIA_CONSTEXPR auto partial_y(const int i, const Point &) ->
            typename Traits<Point>::Scalar {
            switch (i) {
                case 0: {
                    return -1;
                }
                case 1: {
                    return 0;
                }
                case 2: {
                    return 1;
                }
                default: {
                    return 0.0;
                }
            }

            return 0.0;
        }

        template <typename Point>
        UTOPIA_INLINE_FUNCTION static UTOPIA_CONSTEXPR auto partial_x(const int i, const Point &) ->
            typename Traits<Point>::Scalar {
            switch (i) {
                case 0: {
                    return -1.;
                }
                case 1: {
                    return 1;
                }
                case 2: {
                    return 0;
                }
                default: {
                    UTOPIA_DEVICE_ASSERT_CXX14(false);
                    return 0.0;
                }
            }

            return 0.0;
        }

        // space-time mixed derivative
        template <typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION static void grad_x_partial_t(const int /*i*/, const Point &, Deriv &dst) {
            // project t coordinates to 0
            UTOPIA_DEVICE_ASSERT_CXX14(dst.size() == 1);
            dst[0] = 0;
        }

        template <typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION static void grad(const int i, const Point &, Grad &g) {
            switch (i) {
                case 0: {
                    g[0] = -1.;
                    g[1] = -1.;
                    return;
                }
                case 1: {
                    g[0] = 1;
                    g[1] = 0;
                    return;
                }
                case 2: {
                    g[0] = 0;
                    g[1] = 1;
                    return;
                }
                default: {
                    g[0] = 0.0;
                    g[1] = 0.0;
                    return;
                }
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_UNIFORM_TRI3_HPP