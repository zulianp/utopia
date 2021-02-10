#ifndef UTOPIA_SIMD_UNIFORM_HEX8_HPP
#define UTOPIA_SIMD_UNIFORM_HEX8_HPP

#include "utopia_DeviceNumber.hpp"
#include "utopia_Elem.hpp"
#include "utopia_Literal.hpp"
#include "utopia_MemType.hpp"
#include "utopia_Quad4.hpp"
#include "utopia_UniformHex8.hpp"
#include "utopia_Views.hpp"

namespace utopia {

    namespace simd_v2 {

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
            using GradValue = utopia::simd_v2::Vector<Scalar, Dim>;
            using STGradX = utopia::simd_v2::Vector<Scalar, Dim - 1>;
            using FunValue = Scalar;
            // using Side = utopia::Quad4<Scalar, Dim>;

            // template <typename IntArray>
            // UTOPIA_INLINE_FUNCTION void side_idx(const std::size_t &i, IntArray &local_side_idx) const {
            //     switch (i) {
            //         case 0: {
            //             local_side_idx[0] = 0;
            //             local_side_idx[1] = 1;
            //             local_side_idx[2] = 5;
            //             local_side_idx[3] = 4;
            //             return;
            //         }
            //         case 1: {
            //             local_side_idx[0] = 1;
            //             local_side_idx[1] = 2;
            //             local_side_idx[2] = 6;
            //             local_side_idx[3] = 5;
            //             return;
            //         }
            //         case 2: {
            //             local_side_idx[0] = 3;
            //             local_side_idx[1] = 2;
            //             local_side_idx[2] = 6;
            //             local_side_idx[3] = 7;
            //             return;
            //         }
            //         case 3: {
            //             local_side_idx[0] = 0;
            //             local_side_idx[1] = 4;
            //             local_side_idx[2] = 7;
            //             local_side_idx[3] = 3;
            //             return;
            //         }
            //         case 4: {
            //             local_side_idx[0] = 0;
            //             local_side_idx[1] = 1;
            //             local_side_idx[2] = 2;
            //             local_side_idx[3] = 3;
            //             return;
            //         }
            //         case 5: {
            //             local_side_idx[0] = 4;
            //             local_side_idx[1] = 5;
            //             local_side_idx[2] = 6;
            //             local_side_idx[3] = 7;
            //             return;
            //         }

            //         default: {
            //             UTOPIA_DEVICE_ASSERT(false);
            //         }
            //     }
            // }

            // UTOPIA_INLINE_FUNCTION void side(const std::size_t &i, Side &side) const {
            //     switch (i) {
            //         case 0: {
            //             node(0, side.node(0));
            //             node(1, side.node(1));
            //             node(5, side.node(2));
            //             node(4, side.node(3));
            //             break;
            //         }
            //         case 1: {
            //             node(1, side.node(0));
            //             node(2, side.node(1));
            //             node(6, side.node(2));
            //             node(5, side.node(3));
            //             break;
            //         }
            //         case 2: {
            //             node(3, side.node(0));
            //             node(2, side.node(1));
            //             node(6, side.node(2));
            //             node(7, side.node(3));
            //             break;
            //         }
            //         case 3: {
            //             node(0, side.node(0));
            //             node(4, side.node(1));
            //             node(7, side.node(2));
            //             node(3, side.node(3));
            //             break;
            //         }
            //         case 4: {
            //             node(0, side.node(0));
            //             node(1, side.node(1));
            //             node(2, side.node(2));
            //             node(3, side.node(3));
            //             break;
            //         }
            //         case 5: {
            //             node(4, side.node(0));
            //             node(5, side.node(1));
            //             node(6, side.node(2));
            //             node(7, side.node(3));
            //             break;
            //         }

            //         default: {
            //             UTOPIA_DEVICE_ASSERT(false);
            //             break;
            //         }
            //     }

            //     side.init(true);
            // }

            template <typename Point>
            UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> decltype(RefHex8::fun(i, p)) {
                return RefHex8::fun(i, p);
            }

            template <typename Point>
            UTOPIA_INLINE_FUNCTION void node(const std::size_t &i, Point &p) const {
                p.copy(translation_);

                switch (i) {
                    case 0: {
                        return;
                    }
                    case 1: {
                        p.add(0, h_[0]);
                        return;
                    }
                    case 2: {
                        p.add(0, h_[0]);
                        p.add(1, h_[1]);
                        return;
                    }
                    case 3: {
                        p.add(1, h_[1]);
                        return;
                    }
                    case 4: {
                        p.add(2, h_[2]);
                        return;
                    }
                    case 5: {
                        p.add(0, h_[0]);
                        p.add(2, h_[2]);
                        return;
                    }
                    case 6: {
                        p.add(0, h_[0]);
                        p.add(1, h_[1]);
                        p.add(2, h_[2]);
                        return;
                    }
                    case 7: {
                        p.add(1, h_[1]);
                        p.add(2, h_[2]);
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
                g.divide(0, h_[0]);
                g.divide(1, h_[1]);
                g.divide(2, h_[2]);
            }

            // // space-time spatial gradient
            // template <typename Point>
            // UTOPIA_INLINE_FUNCTION auto partial_t(const int i, const Point &p) -> typename Traits<Point>::Scalar {
            //     return RefHex8::partial_t(i, p) / h_[2];
            // }

            // // space-time spatial gradient
            // template <typename Point, typename Deriv>
            // UTOPIA_INLINE_FUNCTION void grad_x(const int i, const Point &p, Deriv &dst) {
            //     RefHex8::grad_x(i, p, dst);
            //     dst[0] /= h_[0];
            //     dst[1] /= h_[1];
            // }

            // /// space-time mixed derivative \nabla_x \partial_t \phi(x, t)
            // template <typename Point, typename Deriv>
            // UTOPIA_INLINE_FUNCTION void grad_x_partial_t(const int i, const Point &p, Deriv &dst) {
            //     RefHex8::grad_x_partial_t(i, p, dst);
            //     dst[0] /= (h_[0] * h_[2]);
            //     dst[1] /= (h_[1] * h_[2]);
            // }

            UTOPIA_INLINE_FUNCTION constexpr static bool is_affine() { return true; }

            UTOPIA_INLINE_FUNCTION constexpr static Scalar reference_measure() { return 1.0; }

            UTOPIA_INLINE_FUNCTION Scalar measure() const {
                UTOPIA_DEVICE_ASSERT(h_[0] * h_[1] * h_[2] > 0.0);
                return h_[0] * h_[1] * h_[2];
            }

            template <typename RefPoint, typename PhysicalPoint>
            UTOPIA_INLINE_FUNCTION void point(const RefPoint &in, PhysicalPoint &out) const {
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

            UTOPIA_INLINE_FUNCTION UniformHex8() { h_.set(0.0); }

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
    }  // namespace simd_v2
}  // namespace utopia

#endif  // UTOPIA_SIMD_UNIFORM_HEX8_HPP
