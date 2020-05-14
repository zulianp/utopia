#ifndef UTOPIA_QUAD_4_HPP
#define UTOPIA_QUAD_4_HPP

#include "utopia_UniformQuad4.hpp"
#include "utopia_DeviceOperations.hpp"

namespace utopia {

    template<class Scalar_, int PhysicalDim>
    class Quad4 {
    public:
        using Scalar = Scalar_;
        static const int Dim        = 1;
        static const int NNodes     = 4;
        static const int NFunctions = 4;
        static const int Order = 1;

        using Point         = utopia::StaticVector<Scalar, PhysicalDim>;
        using RefPoint      = utopia::StaticVector<Scalar, Dim>;
        using GradValue     = utopia::StaticVector<Scalar, PhysicalDim>;
        using STGradX       = utopia::StaticVector<Scalar, Dim-1>;
        using FunValue      = Scalar;
        using MemType       = utopia::Varying<>;
        using Side          = utopia::Edge2<Scalar, PhysicalDim>;
        using Jacobian      = utopia::StaticMatrix<Scalar, PhysicalDim, 2>;
        using InverseJacobian = utopia::StaticMatrix<Scalar, 2, PhysicalDim>;

        //IMPLEMENT ME

        virtual ~Quad4() = default;

        template<typename RefPointT>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const RefPointT &p) -> decltype(RefQuad4::fun(i, p))
        {
          return RefQuad4::fun(i, p);
        }

        template<typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION void grad(const int i, const RefPoint &p, Grad &g) const
        {
            UTOPIA_DEVICE_ASSERT(is_affine_);
            StaticVector<Scalar, 2> g_x;
            RefQuad4::grad(i, p, g_x);
            g = affine_jacobian_inverse_ * g_x;
        }

        template<typename Point>
        UTOPIA_INLINE_FUNCTION void node(const std::size_t &i, Point &p) const
        {
            p.copy(nodes_[i]);
        }

        UTOPIA_INLINE_FUNCTION const Point & node(const std::size_t &i) const
        {
            return nodes_[i];
        }

        UTOPIA_INLINE_FUNCTION Point & node(const std::size_t &i)
        {
            return nodes_[i];
        }

        template<typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void centroid(PhysicalPoint &out) const
        {
            out = 0.25 * (nodes_[0] + nodes_[1] + nodes_[2] + nodes_[3]);
        }

        UTOPIA_INLINE_FUNCTION bool is_affine() const
        {
            return is_affine_;
        }

        UTOPIA_INLINE_FUNCTION constexpr static Scalar reference_measure()
        {
            return 1.0;
        }

        UTOPIA_INLINE_FUNCTION Scalar measure() const
        {
            return measure_;
        }

        // template<typename PhysicalPoint>
        // UTOPIA_INLINE_FUNCTION bool contains(const PhysicalPoint &p, const Scalar &tol = device::epsilon<Scalar>()) const
        // {
        //     Point t = (nodes_[1] - nodes_[0])/h_;
        //     Point n;
        //     n[0] = -t[1];
        //     n[1] =  t[0];

        //     Point v = p - nodes_[0];
        //     if(dot(v, n) > tol*10) {
        //         return false;
        //     }

        //     const Scalar proj = dot(n, t);
        //     if(proj < -tol) {
        //         return false;
        //     }

        //     if(proj > h_ + tol) {
        //         return false;
        //     }

        //     return true;
        // }

        template<typename RefPointT, typename PhysicalPointT>
        UTOPIA_INLINE_FUNCTION void point(const RefPointT &in, PhysicalPointT &out) const
        {
            out = RefQuad4::fun(0, in) * nodes_[0] + RefQuad4::fun(1, in) * nodes_[1] + RefQuad4::fun(2, in) * nodes_[2] + RefQuad4::fun(3, in) * nodes_[3];
        }

        UTOPIA_INLINE_FUNCTION Quad4(const bool &is_affine = true) : measure_(0), is_affine_(is_affine) {}

        void init(const bool is_affine)
        {
            is_affine_ = is_affine;

            if(is_affine_) {
                affine_jacobian_.set_col(0, nodes_[1] - nodes_[0]);
                affine_jacobian_.set_col(1, nodes_[3] - nodes_[0]);

                measure_ = det(affine_jacobian_);

                affine_jacobian_inverse_ = inv(affine_jacobian_);

            } else {
                UTOPIA_DEVICE_ASSERT(false);
            }
        }

        template<typename PhysicalPointT, typename RefPointT>
        UTOPIA_INLINE_FUNCTION void inverse_transform(const PhysicalPointT &in, RefPointT &out) const
        {
            UTOPIA_DEVICE_ASSERT(is_affine_);
            out = affine_jacobian_inverse_ * (in - nodes_[0]);
        }

        UTOPIA_INLINE_FUNCTION constexpr static int n_nodes()
        {
            return NNodes;
        }

        UTOPIA_INLINE_FUNCTION constexpr static int n_functions()
        {
            return NFunctions;
        }

    private:
        Point nodes_[4];
        Scalar measure_;
        Jacobian affine_jacobian_;
        InverseJacobian affine_jacobian_inverse_;
        bool is_affine_;

    };
}

#endif //UTOPIA_QUAD_4_HPP
