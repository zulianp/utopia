#ifndef UTOPIA_REF_EDGE_AND_EDGE_2_HPP
#define UTOPIA_REF_EDGE_AND_EDGE_2_HPP

#include "utopia_Views.hpp"
#include "utopia_DeviceNumber.hpp"
#include "utopia_MemType.hpp"
#include "utopia_Elem.hpp"
#include "utopia_Node1.hpp"


namespace utopia {

    class RefEdge2 {
    public:
        template<typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> typename Traits<Point>::Scalar
        {
            const auto x = p[0];
            switch(i) {
                case 0: { return (1 - x); }
                case 1: { return x; }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return 0.0;
                }
            }
        }

        template<typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION static void grad(const int i, const Point &p, Grad &g)
        {
            const auto x = p[0];
            switch(i) {
                case 0: { g[0] = -1; return; }
                case 1: { g[0] =  1; return; }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                }
            }
        }

        template<typename Point, typename Values>
        UTOPIA_INLINE_FUNCTION static void fun(const Point &p, Values &values)
        {
            const auto x = p[0];
            values[0] = (1 - x);
            values[1] = x;
        }
    };

    template<typename Scalar_, int PhysicalDim>
    class Edge2 : public Elem {
    public:
        using Scalar = Scalar_;
        static const int Dim = PhysicalDim;
        static const int NNodes = 2;
        static const int NSides = 2;
        static const int NFunctions = 2;
        using Point = utopia::StaticVector<Scalar, Dim>;
        using GradValue = utopia::StaticVector<Scalar, Dim>;
        using STGradX   = utopia::StaticVector<Scalar, Dim-1>;
        using FunValue  = Scalar;
        using MemType   = utopia::Varying<>;
        using Side      = utopia::Node1<Scalar, Dim>;

        virtual ~Edge2() {}

        template<typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> decltype(RefEdge2::fun(i, p))
        {
          return RefEdge2::fun(i, p);
        }

        template<typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION void grad(const int i, const Point &p, Grad &g) const
        {
            Scalar g_x[1] = {0.0};
            RefEdge2::grad(i, p, g_x);
            g = (nodes_[1] - nodes_[0]) * (g_x[0] * h_ * h_);
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
            out = 0.5 * (nodes_[0] + nodes_[1]);
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
            return h_;
        }

        template<typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION bool contains(const PhysicalPoint &p, const Scalar &tol = device::epsilon<Scalar>()) const
        {
            Point t = (nodes_[1] - nodes_[0])/h_;
            Point n;
            n[0] = -t[1];
            n[1] =  t[0];

            Point v = p - nodes_[0];
            if(dot(v, n) > tol*10) {
                return false;
            }

            const Scalar proj = dot(n, t);
            if(proj < -tol) {
                return false;
            }

            if(proj > h_ + tol) {
                return false;
            }

            return true;
        }

        template<typename RefPoint, typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void point(const RefPoint &in, PhysicalPoint &out) const
        {
            out = in[0] * nodes_[1] + nodes_[0];
        }

        UTOPIA_INLINE_FUNCTION Edge2() {}

        UTOPIA_INLINE_FUNCTION Edge2(const Point &p1, const Point &p2)
        {
            init(p1, p2);
        }

        UTOPIA_INLINE_FUNCTION void init(const Point &p1, const Point &p2)
        {
            nodes_[0].copy(p1);
            nodes_[1].copy(p2);
            init();
        }

        void init()
        {
            h_ = norm2(nodes_[1] - nodes_[0]);
        }

        template<typename PhysicalPoint, typename RefPoint>
        UTOPIA_INLINE_FUNCTION void inverse_transform(const PhysicalPoint &in, RefPoint &out) const
        {
            out[0] = norm2(in - nodes_[0])/h_;
        }

        UTOPIA_INLINE_FUNCTION constexpr static int n_nodes()
        {
            return NNodes;
        }

        UTOPIA_INLINE_FUNCTION constexpr static int n_functions()
        {
            return NFunctions;
        }

        // UTOPIA_INLINE_FUNCTION const Point &translation() const
        // {
        //     return nodes_[0];
        // }

    private:
        Point nodes_[2];
        Scalar h_;
    };

}

#endif //UTOPIA_REF_EDGE_AND_EDGE_2_HPP
