#ifndef UTOPIA_UNIFORM_QUAD4_HPP
#define UTOPIA_UNIFORM_QUAD4_HPP

#include "utopia_Views.hpp"

namespace utopia {

    class RefQuad4 {
    public:
        template<typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> decltype(p[0])
        {
            const auto x = p[0];
            const auto y = p[1];

            switch(i) {
                case 0: { return (1 - x) * (1 - y); }
                case 1: { return x * (1 - y); }
                case 2: { return x * y; }
                case 3: { return (1 - x) * y; }
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
            const auto y = p[1];

            switch(i)
            {
                case 0:
                {
                    g[0] = y - 1.;
                    g[1] = x - 1.;
                    return;
                }
                case 1:
                {
                    g[0] = 1 - y;
                    g[1] = -x;
                    return;
                }
                case 2:
                {
                    g[0] = y;
                    g[1] = x;
                    return;
                }
                case 3:
                {
                    g[0] = -y;
                    g[1] = (1 - x);
                    return;
                }
                default:
                {
                    g[0] = 0.0;
                    g[1] = 0.0;
                    return;
                }
            }
        }

        template<typename Point, typename Values>
        UTOPIA_INLINE_FUNCTION static void fun(const Point &p, Values &values)
        {
            const auto x = p[0];
            const auto y = p[1];

            values[0] = (1 - x) * (1 - y);
            values[1] = x * (1 - y);
            values[2] = x * y;
            values[3] = (1 - x) * y;
        }

        template<typename Point, typename Values>
        UTOPIA_INLINE_FUNCTION static void grad(const Point &p, Values &values)
        {
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

    template<typename Scalar_>
    class UniformQuad4 final {
    public:
        using Scalar = Scalar_;
        using MemType = Uniform<>;
        // using DiscretizationType = FE;
        static const int Dim = 2;
        static const int NNodes = 4;

        using NodeIndexView = utopia::ArrayView<std::size_t, NNodes>;

        template<typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> decltype(RefQuad4::fun(i, p))
        {
          return RefQuad4::fun(i, p);
        }

        template<typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION void grad(const int i, const Point &p, Grad &g) const
        {
            RefQuad4::grad(i, p, g);
            g[0] /= h_[0];
            g[1] /= h_[1];
        }

        UTOPIA_INLINE_FUNCTION constexpr static bool is_affine()
        {
            return true;
        }

        UTOPIA_INLINE_FUNCTION constexpr static Scalar reference_measure()
        {
            return 1.0;
        }

        // template<typename Point, typename Values>
        // UTOPIA_INLINE_FUNCTION static void fun(const Point &p, Values &values)
        // {
        //     RefQuad4::fun(p, values);
        // }

        // template<typename Point, typename Values>
        // UTOPIA_INLINE_FUNCTION void grad(const Point &p, Values &values) const
        // {
        //     RefQuad4::grad(p, values);
        //     for(int i = 0; i < 4; ++i)
        //     {
        //         values(i, 0) /= h_[0];
        //         values(i, 1) /= h_[1];
        //     }
        // }

        template<typename RefPoint, typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void point(const RefPoint &in, PhysicalPoint &out) const
        {
            out[0] = in[0] * h_[0] + translation_[0];
            out[1] = in[1] * h_[1] + translation_[1];
        }

        UTOPIA_INLINE_FUNCTION UniformQuad4(const Scalar &hx, const Scalar &hy)
        {
            h_[0] = hx;
            h_[1] = hy;
        }

        UTOPIA_INLINE_FUNCTION UniformQuad4()
        {
            h_[0] = 0.0;
            h_[1] = 0.0;
        }

        template<class H>
        UTOPIA_INLINE_FUNCTION void set(
            const StaticVector2<Scalar> &translation,
            const H &h)
        {
            translation_(0) = translation(0);
            translation_(1) = translation(1);

            h_[0] = h[0];
            h_[1] = h[1];
        }

        UTOPIA_INLINE_FUNCTION NodeIndexView &nodes()
        {
            return nodes_;
        }

        UTOPIA_INLINE_FUNCTION const NodeIndexView &nodes() const
        {
            return nodes_;
        }

        UTOPIA_INLINE_FUNCTION const std::size_t &node(const std::size_t &i) const
        {
            return nodes_[i];
        }

        UTOPIA_INLINE_FUNCTION constexpr static int n_nodes()
        {
            return NNodes;
        }

    private:
        Scalar h_[2];
        StaticVector2<Scalar> translation_;
        NodeIndexView nodes_;
    };

}

#endif //UTOPIA_UNIFORM_QUAD4_HPP
