#ifndef UTOPIA_UNIFORM_QUAD4_HPP
#define UTOPIA_UNIFORM_QUAD4_HPP

#include "utopia_DeviceNumber.hpp"
#include "utopia_Edge2.hpp"
#include "utopia_Elem.hpp"
#include "utopia_MemType.hpp"
#include "utopia_Views.hpp"

namespace utopia {

    class RefQuad4 {
    public:
        template <typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> typename Traits<Point>::Scalar {
            const auto x = p[0];
            const auto y = p[1];

            switch (i) {
                case 0: {
                    return (1 - x) * (1 - y);
                }
                case 1: {
                    return x * (1 - y);
                }
                case 2: {
                    return x * y;
                }
                case 3: {
                    return (1 - x) * y;
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
            const auto x = p[0];

            switch (i) {
                case 0: {
                    return x - 1.;
                }
                case 1: {
                    return -x;
                }
                case 2: {
                    return x;
                }
                case 3: {
                    return (1 - x);
                }
                default: {
                    return 0.0;
                }
            }
        }

        // space-time spatial gradient
        template <typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION static void grad_x(const int i, const Point &p, Deriv &dst) {
            UTOPIA_DEVICE_ASSERT(dst.size() == 1);

            // const auto x = p[0];
            const auto t = p[1];

            switch (i) {
                case 0: {
                    dst[0] = t - 1.;
                    return;
                }
                case 1: {
                    dst[0] = 1 - t;
                    return;
                }
                case 2: {
                    dst[0] = t;
                    return;
                }
                case 3: {
                    dst[0] = -t;
                    return;
                }
                default: {
                    dst[0] = 0.0;
                    return;
                }
            }
        }

        // space-time mixed derivative
        template <typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION static void grad_x_partial_t(const int i, const Point & /*p*/, Deriv &dst) {
            // project t coordinates to 0
            UTOPIA_DEVICE_ASSERT(dst.size() == 1);

            switch (i) {
                case 0: {
                    dst[0] = 1.0;
                    return;
                }
                case 1: {
                    dst[0] = -1.0;
                    return;
                }
                case 2: {
                    dst[0] = 1.0;
                    return;
                }
                case 3: {
                    dst[0] = -1.0;
                    return;
                }
                default: {
                    dst[0] = 0.0;
                    return;
                }
            }
        }

        template <typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION static void grad(const int i, const Point &p, Grad &g) {
            const auto x = p[0];
            const auto y = p[1];

            switch (i) {
                case 0: {
                    g[0] = y - 1.;
                    g[1] = x - 1.;
                    return;
                }
                case 1: {
                    g[0] = 1 - y;
                    g[1] = -x;
                    return;
                }
                case 2: {
                    g[0] = y;
                    g[1] = x;
                    return;
                }
                case 3: {
                    g[0] = -y;
                    g[1] = (1 - x);
                    return;
                }
                default: {
                    g[0] = 0.0;
                    g[1] = 0.0;
                    return;
                }
            }
        }

        template <typename Point, typename Values>
        UTOPIA_INLINE_FUNCTION static void fun(const Point &p, Values &values) {
            const auto x = p[0];
            const auto y = p[1];

            values[0] = (1 - x) * (1 - y);
            values[1] = x * (1 - y);
            values[2] = x * y;
            values[3] = (1 - x) * y;
        }

        template <typename Point, typename Values>
        UTOPIA_INLINE_FUNCTION static void grad(const Point &p, Values &values) {
            const auto x = p[0];
            const auto y = p[1];

            values(0, 0) = y - 1.;
            values(0, 1) = x - 1.;

            values(1, 0) = 1 - y;
            values(1, 1) = -x;

            values(2, 0) = y;
            values(2, 1) = x;

            values(3, 0) = -y;
            values(3, 1) = (1 - x);
        }
    };

    template <typename Scalar_>
    class UniformQuad4 : public Elem {
    public:
        using Scalar = Scalar_;
        using MemType = Uniform<>;
        static const int Dim = 2;
        static const int NNodes = 4;
        static const int NSides = 4;
        static const int NFunctions = 4;
        static const int Order = 1;

        using Side = utopia::Edge2<Scalar, Dim>;
        using Point = utopia::StaticVector<Scalar, Dim>;
        using GradValue = utopia::StaticVector<Scalar, Dim>;
        using STGradX = utopia::StaticVector<Scalar, Dim - 1>;
        using FunValue = Scalar;

        template <typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> decltype(RefQuad4::fun(i, p)) {
            return RefQuad4::fun(i, p);
        }

        template <typename Point>
        UTOPIA_INLINE_FUNCTION void node(const std::size_t &i, Point &p) const {
            p[0] = translation_[0];
            p[1] = translation_[1];

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
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return;
                }
            }
        }

        template <typename IntArray>
        UTOPIA_INLINE_FUNCTION void side_idx(const std::size_t &i, IntArray &local_side_idx) const {
            switch (i) {
                case 0: {
                    local_side_idx[0] = 0;
                    local_side_idx[1] = 1;
                    return;
                }
                case 1: {
                    local_side_idx[0] = 1;
                    local_side_idx[1] = 2;
                    return;
                }
                case 2: {
                    local_side_idx[0] = 2;
                    local_side_idx[1] = 3;
                    return;
                }
                case 3: {
                    local_side_idx[0] = 3;
                    local_side_idx[1] = 0;
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
                    break;
                }
                case 1: {
                    node(1, side.node(0));
                    node(2, side.node(1));
                    break;
                }
                case 2: {
                    node(2, side.node(0));
                    node(3, side.node(1));
                    break;
                }
                case 3: {
                    node(3, side.node(0));
                    node(0, side.node(1));
                    break;
                }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    break;
                }
            }

            side.init();
        }

        template <typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void centroid(PhysicalPoint &out) const {
            for (int i = 0; i < Dim; ++i) {
                out[i] = translation_[i] + h_[i] / 2.0;
            }
        }

        // space-time spatial gradient
        template <typename Point>
        UTOPIA_INLINE_FUNCTION auto partial_t(const int i, const Point &p) -> typename Traits<Point>::Scalar {
            return RefQuad4::partial_t(i, p) / h_[1];
        }

        // space-time spatial gradient
        template <typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION void grad_x(const int i, const Point &p, Deriv &dst) {
            RefQuad4::grad_x(i, p, dst);
            dst[0] /= h_[0];
        }

        /// space-time mixed derivative \nabla_x \partial_t \phi(x, t)
        template <typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION void grad_x_partial_t(const int i, const Point &p, Deriv &dst) {
            RefQuad4::grad_x_partial_t(i, p, dst);
            dst[0] /= (h_[0] * h_[1]);
        }

        template <typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION void grad(const int i, const Point &p, Grad &g) const {
            RefQuad4::grad(i, p, g);
            g[0] /= h_[0];
            g[1] /= h_[1];
        }

        UTOPIA_INLINE_FUNCTION constexpr static bool is_affine() { return true; }

        UTOPIA_INLINE_FUNCTION constexpr static Scalar reference_measure() { return 1.0; }

        UTOPIA_INLINE_FUNCTION Scalar measure() const { return h_[0] * h_[1]; }

        template <typename RefPoint, typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void point(const RefPoint &in, PhysicalPoint &out) const {
            out[0] = in[0] * h_[0] + translation_[0];
            out[1] = in[1] * h_[1] + translation_[1];
        }

        template <typename PhysicalPoint, typename RefPoint>
        UTOPIA_INLINE_FUNCTION void inverse_transform(const PhysicalPoint &in, RefPoint &out) const {
            out[0] = (in[0] - translation_[0]) / h_[0];
            out[1] = (in[1] - translation_[1]) / h_[1];
        }

        UTOPIA_INLINE_FUNCTION UniformQuad4(const Scalar &hx, const Scalar &hy) {
            h_[0] = hx;
            h_[1] = hy;
        }

        UTOPIA_INLINE_FUNCTION UniformQuad4() {
            h_[0] = 0.0;
            h_[1] = 0.0;
        }

        template <class H>
        UTOPIA_INLINE_FUNCTION void set(const StaticVector2<Scalar> &translation, const H &h) {
            translation_(0) = translation(0);
            translation_(1) = translation(1);

            h_[0] = h[0];
            h_[1] = h[1];
        }

        bool contains(const Point &p, const Scalar tol = 0.0) const {
            for (int i = 0; i < Dim; ++i) {
                if ((p[i] + tol) <= translation_[i] || (p[i] - tol) > (translation_[i] + h_[i])) return false;
            }

            return true;
        }

        bool intersect_line(const Point &a, const Point &b, Point &a_out, Point &b_out) const {
            // cheap detect
            for (int d = 0; d < Dim; ++d) {
                if (device::max(a[d], b[d]) <= translation_[d] || device::min(a[d], b[d]) >= translation_[d] + h_[d]) {
                    return false;
                }
            }

            // copy content
            a_out.copy(a);
            b_out.copy(b);

            // left-boundary
            Scalar dist_a = a_out[0] - translation_[0];
            Scalar dist_b = b_out[0] - translation_[0];
            Scalar len = b_out[0] - a_out[0];

            bool intersects = false;

            if (device::signbit(dist_a) != device::signbit(dist_b)) {
                intersects = true;

                const Scalar ratio_a = dist_b / len;
                const Scalar ratio_b = -dist_a / len;

                if (dist_a <= 0) {
                    a_out = ratio_a * a_out + ratio_b * b_out;
                } else {
                    b_out = ratio_a * a_out + ratio_b * b_out;
                }

                len = b_out[0] - a_out[0];
            } else {
                if (dist_a <= 0 && dist_b <= 0) {
                    return false;
                }
            }

            // right-boundary
            dist_a = a_out[0] - translation_[0] - h_[0];
            dist_b = b_out[0] - translation_[0] - h_[0];

            if (device::signbit(dist_a) != device::signbit(dist_b)) {
                intersects = true;

                const Scalar ratio_a = dist_b / len;
                const Scalar ratio_b = -dist_a / len;

                if (dist_a >= 0) {
                    a_out = ratio_a * a_out + ratio_b * b_out;
                } else {
                    b_out = ratio_a * a_out + ratio_b * b_out;
                }

                len = b_out[0] - a_out[0];
            } else {
                if (dist_a >= 0 && dist_b >= 0) {
                    return false;
                }
            }

            // bottom-boundary
            dist_a = a_out[1] - translation_[1];
            dist_b = b_out[1] - translation_[1];
            len = b_out[1] - a_out[1];

            if (device::signbit(dist_a) != device::signbit(dist_b)) {
                intersects = true;

                const Scalar ratio_a = dist_b / len;
                const Scalar ratio_b = -dist_a / len;

                if (dist_a <= 0) {
                    a_out = ratio_a * a_out + ratio_b * b_out;
                } else {
                    b_out = ratio_a * a_out + ratio_b * b_out;
                }

                len = b_out[1] - a_out[1];
            } else {
                if (dist_a <= 0 && dist_b <= 0) {
                    return false;
                }
            }

            // top-boundary
            dist_a = a_out[1] - translation_[1] - h_[1];
            dist_b = b_out[1] - translation_[1] - h_[1];

            if (device::signbit(dist_a) != device::signbit(dist_b)) {
                intersects = true;

                const Scalar ratio_a = dist_b / len;
                const Scalar ratio_b = -dist_a / len;

                if (dist_a >= 0) {
                    a_out = ratio_a * a_out + ratio_b * b_out;
                } else {
                    b_out = ratio_a * a_out + ratio_b * b_out;
                }
            } else {
                if (dist_a >= 0 && dist_b >= 0) {
                    return false;
                }
            }

            // disp("---------------");
            // disp("line");
            // disp(a);
            // disp(b);
            // disp("quad");
            // disp(translation_);
            // Point top = translation_ + h_;
            // disp(top);
            // disp("result");
            // disp(a_out);
            // disp(b_out);
            // disp("---------------");

            return intersects;
        }

        UTOPIA_INLINE_FUNCTION constexpr static int n_nodes() { return NNodes; }

        UTOPIA_INLINE_FUNCTION const StaticVector2<Scalar> &translation() const { return translation_; }

    private:
        StaticVector2<Scalar> h_;
        StaticVector2<Scalar> translation_;
    };

}  // namespace utopia

#endif  // UTOPIA_UNIFORM_QUAD4_HPP
