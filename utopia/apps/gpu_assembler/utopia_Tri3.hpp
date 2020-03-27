#ifndef UTOPIA_TRI_3_HPP
#define UTOPIA_TRI_3_HPP

#include "utopia_UniformTri3.hpp"
#include "utopia_ElemTraits.hpp"

namespace utopia {

    template<typename Scalar_, int PhysicalDim = 2>
    class Tri3 : public Elem {
    public:
        using Scalar = Scalar_;
        using MemType = Uniform<>;
        static const int Dim        = 2;
        static const int NNodes     = 3;
        static const int NSides     = 3;
        static const int NFunctions = 3;
        static const int Order = 1;

        using Side             = utopia::Edge2<Scalar, PhysicalDim>;
        using Point            = utopia::StaticVector<Scalar, PhysicalDim>;
        using RefPoint         = utopia::StaticVector<Scalar, Dim>;
        using GradValue        = utopia::StaticVector<Scalar, PhysicalDim>;
        using STGradX          = utopia::StaticVector<Scalar, PhysicalDim-1>;
        using FunValue         = Scalar;

        using Jacobian         = utopia::StaticMatrix<Scalar, PhysicalDim, Dim>;
        using JacobianInverse  = utopia::StaticMatrix<Scalar, Dim, PhysicalDim>;

        template<typename Point>
        UTOPIA_INLINE_FUNCTION static Scalar fun(const int i, const Point &p)
        {
            return RefTri3::fun(i, p);
        }

        template<typename Point>
        UTOPIA_INLINE_FUNCTION void node(const std::size_t &i, Point &p) const
        {
            p.copy(translation_);
            auto &&data = jacobian_.raw_type();

            switch(i) {
                case 0: { return; }
                case 1:
                {
                    for(int d = 0; d < PhysicalDim; ++d)  {
                        p(d) += data[d*Dim];
                    }

                    return;
                }
                case 2:
                {
                    for(int d = 0; d < PhysicalDim; ++d)  {
                        p(d) += data[d*Dim + 1];
                    }

                    return;
                }
                default:
                {
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
                    local_side_idx[0] = 0;
                    local_side_idx[1] = 1;
                    return;
                }
                case 1:
                {
                    local_side_idx[0] = 1;
                    local_side_idx[1] = 2;
                    return;
                }
                case 2:
                {
                    local_side_idx[0] = 2;
                    local_side_idx[1] = 0;
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

        template<typename PhysicalPointT>
        UTOPIA_INLINE_FUNCTION void centroid(PhysicalPointT &out) const
        {
            auto &&data = jacobian_.raw_type();

            for(int i = 0; i < Dim; ++i) {
                const int offset = i*PhysicalDim;
                out[i] = (3.0 * translation_[i] + data[offset] + data[offset + 1])/3.0;
            }
        }

        // //space-time spatial gradient
        template<typename Point>
        UTOPIA_INLINE_FUNCTION auto partial_t(const int i, const Point &p) -> typename Traits<Point>::Scalar
        {
            return jacobian_inverse_(1, 0) * RefTri3::partial_x(i, p) + jacobian_inverse_(1, 1) * RefTri3::partial_y(i, p);
        }

        // //space-time spatial gradient
        template<typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION void grad_x(const int i, const Point &p, Deriv &dst)
        {
           dst[0] = jacobian_inverse_(0, 0) * RefTri3::partial_x(i, p) + jacobian_inverse_(0, 1) *  RefTri3::partial_y(i, p);
        }

        // ///space-time mixed derivative \nabla_x \partial_t \phi(x, t)
        template<typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION void grad_x_partial_t(const int , const Point &, Deriv &dst)
        {
            dst[0] = 0;
        }

        template<typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION void grad(const int i, const Point &p, Grad &g) const
        {
            RefPoint g_ref;
            RefTri3::grad(i, p, g_ref);
            g = jacobian_inverse_ * g_ref;
        }

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
            return measure_;
        }

        template<typename RefPointT, typename PhysicalPointT>
        UTOPIA_INLINE_FUNCTION void point(const RefPointT &in, PhysicalPointT &out) const
        {
            out = jacobian_ * in + translation_;
        }

        template<typename PhysicalPointT, typename RefPointT>
        UTOPIA_INLINE_FUNCTION void inverse_transform(const PhysicalPointT &in, RefPointT &out) const
        {
            out = jacobian_inverse_ * (in - translation_);
        }

        UTOPIA_INLINE_FUNCTION Tri3() {}

        template<class P0, class P1, class P2>
        UTOPIA_INLINE_FUNCTION void set(
            const P0 &p0,
            const P1 &p1,
            const P2 &p2
        )
        {
            translation_.copy(p0);

            jacobian_.set_col(0, p1 - translation_);
            jacobian_.set_col(1, p2 - translation_);

            jacobian_inverse_ = inv(jacobian_);

            measure_ = det(jacobian_);
            UTOPIA_DEVICE_ASSERT(measure_>0);
        }

        UTOPIA_INLINE_FUNCTION constexpr static int n_nodes()
        {
            return NNodes;
        }

        UTOPIA_INLINE_FUNCTION const Point &translation() const
        {
            return translation_;
        }

    private:
        Point translation_;

        //FIXME
        #pragma GCC diagnostic warning "-Wuninitialized"
        Jacobian jacobian_;

        #pragma GCC diagnostic warning "-Wuninitialized"
        JacobianInverse jacobian_inverse_;
        Scalar measure_;
    };

    template<typename Scalar, int PhysicalDim>
    struct is_simplex<Tri3<Scalar, PhysicalDim>> {
        static const bool value = true;
    };

}

#endif //UTOPIA_TRI_3_HPP
