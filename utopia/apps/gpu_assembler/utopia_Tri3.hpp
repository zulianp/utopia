#ifndef UTOPIA_TRI_3_HPP
#define UTOPIA_TRI_3_HPP

#include "utopia_UniformTri3.hpp"

namespace utopia {
    template<typename Scalar_>
    class Tri3 : public Elem {
    public:
        using Scalar = Scalar_;
        using MemType = Uniform<>;
        static const int Dim        = 2;
        static const int NNodes     = 3;
        static const int NSides     = 3;
        static const int NFunctions = 3;

        using Side      = utopia::Edge2<Scalar, Dim>;
        using Point     = utopia::StaticVector<Scalar, Dim>;
        using GradValue = utopia::StaticVector<Scalar, Dim>;
        using STGradX   = utopia::StaticVector<Scalar, Dim-1>;
        using FunValue  = Scalar;


        template<typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> decltype(RefTri3::fun(i, p))
        {
          return RefTri3::fun(i, p);
        }

        template<typename Point>
        UTOPIA_INLINE_FUNCTION void node(const std::size_t &i, Point &p) const
        {
            p[0] = translation_[0];
            p[1] = translation_[1];

            switch(i) {
                case 0: { return; }
                case 1: { p += u_; return; }
                case 2: { p += v_; return; }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return;
                }
            }
        }

        template<typename IntArray>
        UTOPIA_INLINE_FUNCTION void side_idx(const std::size_t &i, IntArray &local_side_idx) const
        {
            switch(i) {
                case 0:
                {
                    local_side_idx[0]= 0;
                    local_side_idx[1]= 1;
                    return;
                }
                case 1:
                {
                    local_side_idx[0]= 1;
                    local_side_idx[1]= 2;
                    return;
                }
                case 2:
                {
                    local_side_idx[0]= 2;
                    local_side_idx[1]= 0;
                    return;
                }
                default:
                {
                    UTOPIA_DEVICE_ASSERT(false);
                }
            }
        }


        UTOPIA_INLINE_FUNCTION void side(const std::size_t &i, Side &side) const
        {
            switch(i) {
                case 0:
                {
                    node(0, side.node(0));
                    node(1, side.node(1));
                    break;
                }
                case 1:
                {
                    node(1, side.node(0));
                    node(2, side.node(1));
                    break;
                }
                case 2:
                {
                    node(2, side.node(0));
                    node(0, side.node(1));
                    break;
                }
                default:
                {
                    UTOPIA_DEVICE_ASSERT(false);
                    break;
                }
            }

            side.init();
        }

        template<typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void centroid(PhysicalPoint &out) const
        {
            for(int i = 0; i < Dim; ++i) {
                out[i] = (3.0 * translation_[i] + u_[i] + v_[i])/3.0;
            }
        }


        // //space-time spatial gradient
        // template<typename Point>
        // UTOPIA_INLINE_FUNCTION auto partial_t(const int i, const Point &p) -> typename Traits<Point>::Scalar
        // {
        //     return RefTri3::partial_t(i, p) / h_[1];
        // }

        // //space-time spatial gradient
        // template<typename Point, typename Deriv>
        // UTOPIA_INLINE_FUNCTION void grad_x(const int i, const Point &p, Deriv &dst)
        // {
        //     RefTri3::grad_x(i, p, dst);
        //     dst[0] /= h_[0];
        // }

        // ///space-time mixed derivative \nabla_x \partial_t \phi(x, t)
        // template<typename Point, typename Deriv>
        // UTOPIA_INLINE_FUNCTION void grad_x_partial_t(const int i, const Point &p, Deriv &dst)
        // {
        //     RefTri3::grad_x_partial_t(i, p, dst);
        //     dst[0] /= (h_[0]*h_[1]);
        // }

        // template<typename Point, typename Grad>
        // UTOPIA_INLINE_FUNCTION void grad(const int i, const Point &p, Grad &g) const
        // {
        //     RefTri3::grad(i, p, g);
        //     g[0] /= h_[0];
        //     g[1] /= h_[1];
        // }

        UTOPIA_INLINE_FUNCTION constexpr static bool is_affine()
        {
            return true;
        }

        UTOPIA_INLINE_FUNCTION constexpr static Scalar reference_measure()
        {
            return 1.0;
        }

        UTOPIA_INLINE_FUNCTION Scalar measure() const
        {
            return (u_[0] * v_[1]) - ( u_[1] * v_[0] );
        }

        template<typename RefPoint, typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void point(const RefPoint &in, PhysicalPoint &out) const
        {
            out[0] = in[0] * u_[0] + in[1] * v_[0] + translation_[0];
            out[1] = in[0] * v_[1] + in[1] * v_[1] + translation_[1];
        }

        template<typename PhysicalPoint, typename RefPoint>
        UTOPIA_INLINE_FUNCTION void inverse_transform(const PhysicalPoint &in, RefPoint &out) const
        {
            out[0] = (in[0] - translation_[0])/h_[0];
            out[1] = (in[1] - translation_[1])/h_[1];
        }

        UTOPIA_INLINE_FUNCTION Tri3() {}

        template<class U, class V>
        UTOPIA_INLINE_FUNCTION void set(
            const StaticVector2<Scalar> &translation,
            const U &u,
            const V &v
            )
        {
            translation_(0) = translation(0);
            translation_(1) = translation(1);

            u_ = u;
            v_ = v;
        }

        UTOPIA_INLINE_FUNCTION constexpr static int n_nodes()
        {
            return NNodes;
        }

        UTOPIA_INLINE_FUNCTION const StaticVector2<Scalar> &translation() const
        {
            return translation_;
        }

    private:
        //u = p2 - p1; v = p3 - p1
        StaticVector2<Scalar> u_, v_;
        StaticVector2<Scalar> translation_;
    };
}

#endif //UTOPIA_TRI_3_HPP
