#ifndef UTOPIA_UNIFORM_HEX8_HPP
#define UTOPIA_UNIFORM_HEX8_HPP

#include "utopia_DeviceNumber.hpp"
#include "utopia_Elem.hpp"
#include "utopia_Literal.hpp"
#include "utopia_MemType.hpp"
#include "utopia_Quad4.hpp"
#include "utopia_Views.hpp"

namespace utopia {

    class RefHex8 {
    public:
        template <typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> typename Traits<Point>::Scalar {
            using Scalar = typename Traits<Point>::Scalar;

            const auto x = p[0];
            const auto y = p[1];
            const auto z = p[2];

            const auto one = One<Scalar>::value();

            switch (i) {
                case 0: {
                    return (one - x) * (one - y) * (one - z);
                }
                case 1: {
                    return x * (one - y) * (one - z);
                }
                case 2: {
                    return x * y * (one - z);
                }
                case 3: {
                    return (one - x) * y * (one - z);
                }
                case 4: {
                    return (one - x) * (one - y) * z;
                }
                case 5: {
                    return x * (one - y) * z;
                }
                case 6: {
                    return x * y * z;
                }
                case 7: {
                    return (one - x) * y * z;
                }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return 0.0;
                }
            }
        }

        // space-time spatial gradient
        template <typename Point>
        UTOPIA_INLINE_FUNCTION static auto partial_t(const int i, const Point &p) -> typename Traits<Point>::Scalar {
            using Scalar = typename Traits<Point>::Scalar;

            const Scalar x = p[0];
            const Scalar y = p[1];
            // const Scalar t = p[2];

            const auto one = One<Scalar>::value();

            switch (i) {
                case 0: {
                    return -(one - x) * (one - y);
                }
                case 1: {
                    return -x * (one - y);
                }
                case 2: {
                    return -x * y;
                }
                case 3: {
                    return -(one - x) * y;
                }
                case 4: {
                    return (one - x) * (one - y);
                }
                case 5: {
                    return x * (one - y);
                }
                case 6: {
                    return x * y;
                }
                case 7: {
                    return (one - x) * y;
                }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return 0.0;
                }
            }
        }

        // space-time spatial gradient
        template <typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION static void grad_x(const int i, const Point &p, Deriv &dst) {
            UTOPIA_DEVICE_ASSERT(dst.size() == 2);

            using Scalar = typename Traits<Point>::Scalar;

            const Scalar x = p[0];
            const Scalar y = p[1];
            const Scalar t = p[2];

            switch (i) {
                // f = (1.0 - x) * (1.0 - y) * (1.0 - t);
                case 0: {
                    dst[0] = -(1.0 - y) * (1.0 - t);
                    dst[1] = -(1.0 - x) * (1.0 - t);
                    return;
                }

                // f = x * (1.0 - y) * (1.0 - t);
                case 1: {
                    dst[0] = (1.0 - y) * (1.0 - t);
                    dst[1] = -x * (1.0 - t);
                    return;
                }

                // f = x * y * (1.0 - t);
                case 2: {
                    dst[0] = y * (1.0 - t);
                    dst[1] = x * (1.0 - t);
                    return;
                }

                // f = (1.0 - x) * y * (1.0 - t);
                case 3: {
                    dst[0] = -y * (1.0 - t);
                    dst[1] = (1.0 - x) * (1.0 - t);
                    return;
                }

                // f = (1.0 - x) * (1.0 - y) * t;
                case 4: {
                    dst[0] = -(1.0 - y) * t;
                    dst[1] = -(1.0 - x) * t;
                    return;
                }

                // f = x * (1.0 - y) * t;
                case 5: {
                    dst[0] = (1.0 - y) * t;
                    dst[1] = -x * t;
                    return;
                }

                // f = x * y * t;
                case 6: {
                    dst[0] = y * t;
                    dst[1] = x * t;
                    return;
                }

                // f = (1.0 - x) * y * t;
                case 7: {
                    dst[0] = -y * t;
                    dst[1] = (1.0 - x) * t;
                    return;
                }

                default: {
                    dst[0] = 0.0;
                    dst[1] = 0.0;
                    return;
                }
            }
        }

        // space-time mixed derivative
        template <typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION static void grad_x_partial_t(const int i, const Point &p, Deriv &dst) {
            // project t coordinates to 0
            UTOPIA_DEVICE_ASSERT(dst.size() == 2);

            using Scalar = typename Traits<Point>::Scalar;

            const Scalar x = p[0];
            const Scalar y = p[1];
            // const Scalar t = p[2];

            switch (i) {
                // f = (1.0 - x) * (1.0 - y) * (1.0 - t);
                // f_t = -(1.0 - x) * (1.0 - y)
                case 0: {
                    dst[0] = (1.0 - y);
                    dst[1] = (1.0 - x);
                    return;
                }

                // f = x * (1.0 - y) * (1.0 - t);
                // f_t = -x * (1.0 - y)
                case 1: {
                    dst[0] = -(1.0 - y);
                    dst[1] = x;
                    return;
                }

                // f = x * y * (1.0 - t);
                // f_t = -x * y
                case 2: {
                    dst[0] = -y;
                    dst[1] = -x;
                    return;
                }

                // f = (1.0 - x) * y * (1.0 - t);
                // f_t = -(1.0 - x) * y
                case 3: {
                    dst[0] = y;
                    dst[1] = -(1.0 - x);
                    return;
                }

                // f = (1.0 - x) * (1.0 - y) * t;
                // f_t = (1.0 - x) * (1.0 - y)
                case 4: {
                    dst[0] = -(1.0 - y);
                    dst[1] = -(1.0 - x);
                    return;
                }

                // f = x * (1.0 - y) * t;
                // f_t = x * (1.0 - y)
                case 5: {
                    dst[0] = (1.0 - y);
                    dst[1] = -x;
                    return;
                }

                // f = x * y * t;
                // f_t = x * y
                case 6: {
                    dst[0] = y;
                    dst[1] = x;
                    return;
                }

                // f = (1.0 - x) * y * t;
                // f_t = (1.0 - x) * y
                case 7: {
                    dst[0] = -y;
                    dst[1] = (1.0 - x);
                    return;
                }

                default: {
                    dst[0] = 0.0;
                    dst[1] = 0.0;
                    return;
                }
            }
        }

        template <typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION static void grad(const int i, const Point &p, Grad &g) {
            using Scalar = typename Traits<Point>::Scalar;

            const auto one = One<Scalar>::value();
            const auto zero = Zero<Scalar>::value();

            const Scalar x = p[0];
            const Scalar y = p[1];
            const Scalar z = p[2];

            switch (i) {
                // f = (1.0 - x) * (1.0 - y) * (1.0 - z);
                case 0: {
                    g.set(0, -(one - y) * (one - z));
                    g.set(1, -(one - x) * (one - z));
                    g.set(2, -(one - x) * (one - y));
                    return;
                }

                // f = x * (one - y) * (one - z);
                case 1: {
                    g.set(0, (one - y) * (one - z));
                    g.set(1, -x * (one - z));
                    g.set(2, -x * (one - y));
                    return;
                }

                // f = x * y * (one - z);
                case 2: {
                    g.set(0, y * (one - z));
                    g.set(1, x * (one - z));
                    g.set(2, -x * y);
                    return;
                }

                // f = (one - x) * y * (one - z);
                case 3: {
                    g.set(0, -y * (one - z));
                    g.set(1, (one - x) * (one - z));
                    g.set(2, -(one - x) * y);
                    return;
                }

                // f = (one - x) * (one - y) * z;
                case 4: {
                    g.set(0, -(one - y) * z);
                    g.set(1, -(one - x) * z);
                    g.set(2, (one - x) * (one - y));
                    return;
                }

                // f = x * (one - y) * z;
                case 5: {
                    g.set(0, (one - y) * z);
                    g.set(1, -x * z);
                    g.set(2, x * (one - y));
                    return;
                }

                // f = x * y * z;
                case 6: {
                    g.set(0, y * z);
                    g.set(1, x * z);
                    g.set(2, x * y);
                    return;
                }

                // f = (one - x) * y * z;
                case 7: {
                    g.set(0, -y * z);
                    g.set(1, (one - x) * z);
                    g.set(2, (one - x) * y);
                    return;
                }

                default: {
                    g.set(0, zero);
                    g.set(1, zero);
                    g.set(2, zero);
                    return;
                }
            }
        }
    };

    /*
        ExodusII format:

               7 ----------- 6
              /|            /|
             / |           / |
            4 ----------- 5  |
            |  |          |  |
            |  |          |  |
            |  3 ---------|- 2
            | /           | /
            |/            |/
            0 ----------- 1

    */
    template <typename Scalar_>
    class UniformHex8 : public Elem {
    public:
        using Scalar = Scalar_;
        using MemType = Uniform<>;
        static const int Dim = 3;
        static const int NNodes = 8;
        static const int NSides = 6;
        static const int NFunctions = NNodes;
        static const int Order = 1;

        using Point = utopia::StaticVector<Scalar, Dim>;
        using GradValue = utopia::StaticVector<Scalar, Dim>;
        using STGradX = utopia::StaticVector<Scalar, Dim - 1>;
        using FunValue = Scalar;
        using Side = utopia::Quad4<Scalar, Dim>;

        template <typename IntArray>
        UTOPIA_INLINE_FUNCTION void side_idx(const std::size_t &i, IntArray &local_side_idx) const {
            switch (i) {
                case 0: {
                    local_side_idx[0] = 0;
                    local_side_idx[1] = 1;
                    local_side_idx[2] = 5;
                    local_side_idx[3] = 4;
                    return;
                }
                case 1: {
                    local_side_idx[0] = 1;
                    local_side_idx[1] = 2;
                    local_side_idx[2] = 6;
                    local_side_idx[3] = 5;
                    return;
                }
                case 2: {
                    local_side_idx[0] = 3;
                    local_side_idx[1] = 2;
                    local_side_idx[2] = 6;
                    local_side_idx[3] = 7;
                    return;
                }
                case 3: {
                    local_side_idx[0] = 0;
                    local_side_idx[1] = 4;
                    local_side_idx[2] = 7;
                    local_side_idx[3] = 3;
                    return;
                }
                case 4: {
                    local_side_idx[0] = 0;
                    local_side_idx[1] = 1;
                    local_side_idx[2] = 2;
                    local_side_idx[3] = 3;
                    return;
                }
                case 5: {
                    local_side_idx[0] = 4;
                    local_side_idx[1] = 5;
                    local_side_idx[2] = 6;
                    local_side_idx[3] = 7;
                    return;
                }

                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                }
            }
        }

        UTOPIA_INLINE_FUNCTION void side(const std::size_t &i, Side &side) const {
            switch (i) {
                case 0: {
                    node(0, side.node(0));
                    node(1, side.node(1));
                    node(5, side.node(2));
                    node(4, side.node(3));
                    break;
                }
                case 1: {
                    node(1, side.node(0));
                    node(2, side.node(1));
                    node(6, side.node(2));
                    node(5, side.node(3));
                    break;
                }
                case 2: {
                    node(3, side.node(0));
                    node(2, side.node(1));
                    node(6, side.node(2));
                    node(7, side.node(3));
                    break;
                }
                case 3: {
                    node(0, side.node(0));
                    node(4, side.node(1));
                    node(7, side.node(2));
                    node(3, side.node(3));
                    break;
                }
                case 4: {
                    node(0, side.node(0));
                    node(1, side.node(1));
                    node(2, side.node(2));
                    node(3, side.node(3));
                    break;
                }
                case 5: {
                    node(4, side.node(0));
                    node(5, side.node(1));
                    node(6, side.node(2));
                    node(7, side.node(3));
                    break;
                }

                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    break;
                }
            }

            side.init(true);
        }

        template <typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> decltype(RefHex8::fun(i, p)) {
            return RefHex8::fun(i, p);
        }

        template <typename Point>
        UTOPIA_INLINE_FUNCTION void node(const std::size_t &i, Point &p) const {
            p[0] = translation_[0];
            p[1] = translation_[1];
            p[2] = translation_[2];

            switch (i) {
                case 0: {
                    return;
                }
                case 1: {
                    p[0] += h_[0];
                    return;
                }
                case 2: {
                    p[0] += h_[0];
                    p[1] += h_[1];
                    return;
                }
                case 3: {
                    p[1] += h_[1];
                    return;
                }
                case 4: {
                    p[2] += h_[2];
                    return;
                }
                case 5: {
                    p[2] += h_[2];
                    p[0] += h_[0];
                    return;
                }
                case 6: {
                    p[2] += h_[2];
                    p[0] += h_[0];
                    p[1] += h_[1];
                    return;
                }
                case 7: {
                    p[2] += h_[2];
                    p[1] += h_[1];
                    return;
                }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return;
                }
            }
        }

        template <typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION void grad(const int i, const Point &p, Grad &g) const {
            RefHex8::grad(i, p, g);
            // g[0] /= h_[0];
            // g[1] /= h_[1];
            // g[2] /= h_[2];

            g.divide(0, h_[0]);
            g.divide(1, h_[1]);
            g.divide(2, h_[2]);
        }

        // space-time spatial gradient
        template <typename Point>
        UTOPIA_INLINE_FUNCTION auto partial_t(const int i, const Point &p) -> typename Traits<Point>::Scalar {
            return RefHex8::partial_t(i, p) / h_[2];
        }

        // space-time spatial gradient
        template <typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION void grad_x(const int i, const Point &p, Deriv &dst) {
            RefHex8::grad_x(i, p, dst);
            dst[0] /= h_[0];
            dst[1] /= h_[1];
        }

        /// space-time mixed derivative \nabla_x \partial_t \phi(x, t)
        template <typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION void grad_x_partial_t(const int i, const Point &p, Deriv &dst) {
            RefHex8::grad_x_partial_t(i, p, dst);
            dst[0] /= (h_[0] * h_[2]);
            dst[1] /= (h_[1] * h_[2]);
        }

        UTOPIA_INLINE_FUNCTION constexpr static bool is_affine() { return true; }

        UTOPIA_INLINE_FUNCTION constexpr static Scalar reference_measure() { return 1.0; }

        UTOPIA_INLINE_FUNCTION Scalar measure() const {
            UTOPIA_DEVICE_ASSERT(h_[0] * h_[1] * h_[2] > 0.0);
            return h_[0] * h_[1] * h_[2];
        }

        template <typename RefPoint, typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void point(const RefPoint &in, PhysicalPoint &out) const {
            // out[0] = in[0] * h_[0] + translation_[0];
            // out[1] = in[1] * h_[1] + translation_[1];
            // out[2] = in[2] * h_[2] + translation_[2];

            for (int i = 0; i < 3; ++i) {
                out.set(i, in[i] * h_[i] + translation_[i]);
            }
        }

        template <typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void centroid(PhysicalPoint &out) const {
            for (int i = 0; i < Dim; ++i) {
                out[i] = translation_[i] + h_[i] / 2.0;
            }
        }

        UTOPIA_INLINE_FUNCTION UniformHex8() {
            h_[0] = 0.0;
            h_[1] = 0.0;
            h_[2] = 0.0;
        }

        template <class Tr, class H>
        UTOPIA_INLINE_FUNCTION void set(const Tr &translation, const H &h) {
            translation_[0] = translation[0];
            translation_[1] = translation[1];
            translation_[2] = translation[2];

            h_[0] = h[0];
            h_[1] = h[1];
            h_[2] = h[2];
        }

        UTOPIA_INLINE_FUNCTION constexpr static int n_nodes() { return NNodes; }

        UTOPIA_INLINE_FUNCTION const Point &translation() const { return translation_; }

    private:
        Point h_;
        Point translation_;
    };

}  // namespace utopia

#endif  // UTOPIA_UNIFORM_HEX8_HPP
