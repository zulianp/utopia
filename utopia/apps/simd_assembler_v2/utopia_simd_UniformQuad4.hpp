#ifndef UTOPIA_SIMD_UNIFORM_QUAD4_HPP
#define UTOPIA_SIMD_UNIFORM_QUAD4_HPP

#include "utopia_DeviceNumber.hpp"
#include "utopia_Edge2.hpp"
#include "utopia_Elem.hpp"
#include "utopia_MemType.hpp"
#include "utopia_Views.hpp"

#include "utopia_UniformQuad4.hpp"

namespace utopia {

    namespace simd_v2 {
        template <typename T>
        using UniformQuad4 = utopia::UniformQuad4<T>;

        // template <typename Scalar_>
        // class UniformQuad4 : public Elem {
        // public:
        //     using Scalar = Scalar_;
        //     using MemType = Uniform<>;
        //     static const int Dim = 2;
        //     static const int NNodes = 4;
        //     static const int NSides = 4;
        //     static const int NFunctions = 4;
        //     static const int Order = 1;

        //     using Side = utopia::Edge2<Scalar, Dim>;
        //     using Point = utopia::StaticVector<Scalar, Dim>;
        //     using GradValue = utopia::StaticVector<Scalar, Dim>;
        //     using STGradX = utopia::StaticVector<Scalar, Dim - 1>;
        //     using FunValue = Scalar;

        //     template <typename Point>
        //     UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> decltype(RefQuad4::fun(i, p)) {
        //         return RefQuad4::fun(i, p);
        //     }

        //     template <typename Point>
        //     UTOPIA_INLINE_FUNCTION void node(const std::size_t &i, Point &p) const {
        //         p[0] = translation_[0];
        //         p[1] = translation_[1];

        //         switch (i) {
        //             case 0: {
        //                 return;
        //             }
        //             case 1: {
        //                 p[0] += h_[0];
        //                 return;
        //             }
        //             case 2: {
        //                 p[0] += h_[0];
        //                 p[1] += h_[1];
        //                 return;
        //             }
        //             case 3: {
        //                 p[1] += h_[1];
        //                 return;
        //             }
        //             default: {
        //                 UTOPIA_DEVICE_ASSERT(false);
        //                 return;
        //             }
        //         }
        //     }

        //     template <typename IntArray>
        //     UTOPIA_INLINE_FUNCTION void side_idx(const std::size_t &i, IntArray &local_side_idx) const {
        //         switch (i) {
        //             case 0: {
        //                 local_side_idx[0] = 0;
        //                 local_side_idx[1] = 1;
        //                 return;
        //             }
        //             case 1: {
        //                 local_side_idx[0] = 1;
        //                 local_side_idx[1] = 2;
        //                 return;
        //             }
        //             case 2: {
        //                 local_side_idx[0] = 2;
        //                 local_side_idx[1] = 3;
        //                 return;
        //             }
        //             case 3: {
        //                 local_side_idx[0] = 3;
        //                 local_side_idx[1] = 0;
        //                 return;
        //             }
        //             default: {
        //                 UTOPIA_DEVICE_ASSERT(false);
        //             }
        //         }
        //     }

        //     UTOPIA_INLINE_FUNCTION void side(const std::size_t &i, Side &side) const {
        //         switch (i) {
        //             case 0: {
        //                 node(0, side.node(0));
        //                 node(1, side.node(1));
        //                 break;
        //             }
        //             case 1: {
        //                 node(1, side.node(0));
        //                 node(2, side.node(1));
        //                 break;
        //             }
        //             case 2: {
        //                 node(2, side.node(0));
        //                 node(3, side.node(1));
        //                 break;
        //             }
        //             case 3: {
        //                 node(3, side.node(0));
        //                 node(0, side.node(1));
        //                 break;
        //             }
        //             default: {
        //                 UTOPIA_DEVICE_ASSERT(false);
        //                 break;
        //             }
        //         }

        //         side.init();
        //     }

        //     template <typename PhysicalPoint>
        //     UTOPIA_INLINE_FUNCTION void centroid(PhysicalPoint &out) const {
        //         for (int i = 0; i < Dim; ++i) {
        //             out[i] = translation_[i] + h_[i] / 2.0;
        //         }
        //     }

        //     // space-time spatial gradient
        //     // template <typename Point>
        //     // UTOPIA_INLINE_FUNCTION auto partial_t(const int i, const Point &p) -> typename Traits<Point>::Scalar {
        //     //     return RefQuad4::partial_t(i, p) / h_[1];
        //     // }

        //     // // space-time spatial gradient
        //     // template <typename Point, typename Deriv>
        //     // UTOPIA_INLINE_FUNCTION void grad_x(const int i, const Point &p, Deriv &dst) {
        //     //     RefQuad4::grad_x(i, p, dst);
        //     //     dst[0] /= h_[0];
        //     // }

        //     // /// space-time mixed derivative \nabla_x \partial_t \phi(x, t)
        //     // template <typename Point, typename Deriv>
        //     // UTOPIA_INLINE_FUNCTION void grad_x_partial_t(const int i, const Point &p, Deriv &dst) {
        //     //     RefQuad4::grad_x_partial_t(i, p, dst);
        //     //     dst[0] /= (h_[0] * h_[1]);
        //     // }

        //     template <typename Point, typename Grad>
        //     UTOPIA_INLINE_FUNCTION void grad(const int i, const Point &p, Grad &g) const {
        //         RefQuad4::grad(i, p, g);
        //         g.divide(0, h_[0]);
        //         g.divide(1, h_[1]);
        //     }

        //     UTOPIA_INLINE_FUNCTION constexpr static bool is_affine() { return true; }

        //     UTOPIA_INLINE_FUNCTION constexpr static Scalar reference_measure() { return 1.0; }

        //     UTOPIA_INLINE_FUNCTION Scalar measure() const { return h_[0] * h_[1]; }

        //     template <typename RefPoint, typename PhysicalPoint>
        //     UTOPIA_INLINE_FUNCTION void point(const RefPoint &in, PhysicalPoint &out) const {
        //         for (int i = 0; i < 2; ++i) {
        //             out.set(i, in[i] * h_[i] + translation_[i]);
        //         }
        //     }

        //     template <typename PhysicalPoint, typename RefPoint>
        //     UTOPIA_INLINE_FUNCTION void inverse_transform(const PhysicalPoint &in, RefPoint &out) const {
        //         for (int i = 0; i < 2; ++i) {
        //             out.set(i, (in[i] - translation_[i]) / h_[i]);
        //         }
        //     }

        //     UTOPIA_INLINE_FUNCTION UniformQuad4(const Scalar &hx, const Scalar &hy) {
        //         h_[0] = hx;
        //         h_[1] = hy;
        //     }

        //     UTOPIA_INLINE_FUNCTION UniformQuad4() {
        //         h_[0] = 0.0;
        //         h_[1] = 0.0;
        //     }

        //     template <class Tr, class H>
        //     UTOPIA_INLINE_FUNCTION void set(const Tr &translation, const H &h) {
        //         translation_[0] = translation[0];
        //         translation_[1] = translation[1];

        //         h_[0] = h[0];
        //         h_[1] = h[1];
        //     }

        //     UTOPIA_INLINE_FUNCTION constexpr static int n_nodes() { return NNodes; }

        //     UTOPIA_INLINE_FUNCTION const StaticVector2<Scalar> &translation() const { return translation_; }

        // private:
        //     StaticVector2<Scalar> h_;
        //     StaticVector2<Scalar> translation_;
        // };

    }  // namespace simd_v2
}  // namespace utopia

#endif  // UTOPIA_SIMD_UNIFORM_QUAD4_HPP
