#ifndef UTOPIA_UNIFORM_HEX8_HPP
#define UTOPIA_UNIFORM_HEX8_HPP


#include "utopia_Views.hpp"
#include "utopia_DeviceNumber.hpp"
#include "utopia_MemType.hpp"
#include "utopia_Elem.hpp"

namespace utopia {

    class RefHex8 {
    public:
        template<typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> typename Traits<Point>::Scalar
        {
            using Scalar = typename Traits<Point>::Scalar;

            const Scalar x = p[0];
            const Scalar y = p[1];
            const Scalar z = p[2];

            switch(i) {
                case 0: { return (1.0 - x) * (1.0 - y) * (1.0 - z); }
                case 1: { return x * (1.0 - y) * (1.0 - z); }
                case 2: { return x * y * (1.0 - z); }
                case 3: { return (1.0 - x) * y * (1.0 - z); }
                case 4: { return (1.0 - x) * (1.0 - y) * z; }
                case 5: { return x * (1.0 - y) * z; }
                case 6: { return x * y * z; }
                case 7: { return (1.0 - x) * y * z; }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return 0.0;
                }
            }
        }


        //space-time spatial gradient
        template<typename Point>
        UTOPIA_INLINE_FUNCTION static auto partial_t(const int i, const Point &p) -> typename Traits<Point>::Scalar
        {
            using Scalar = typename Traits<Point>::Scalar;

            const Scalar x = p[0];
            const Scalar y = p[1];
            const Scalar t = p[2];

            switch(i) {
                case 0: { return -(1.0 - x) * (1.0 - y); }
                case 1: { return -x * (1.0 - y); }
                case 2: { return -x * y; }
                case 3: { return -(1.0 - x) * y; }
                case 4: { return (1.0 - x) * (1.0 - y); }
                case 5: { return x * (1.0 - y); }
                case 6: { return x * y; }
                case 7: { return (1.0 - x) * y; }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return 0.0;
                }
            }
        }

        //space-time spatial gradient
        template<typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION static void grad_x(const int i, const Point &p, Deriv &dst)
        {
            UTOPIA_DEVICE_ASSERT(dst.size() == 2);

            using Scalar = typename Traits<Point>::Scalar;

            const Scalar x = p[0];
            const Scalar y = p[1];
            const Scalar t = p[2];

            switch(i)
            {
                // f = (1.0 - x) * (1.0 - y) * (1.0 - t);
                case 0:
                {
                    dst[0] = -(1.0 - y) * (1.0 - t);
                    dst[1] = -(1.0 - x) * (1.0 - t);
                    return;
                }

                // f = x * (1.0 - y) * (1.0 - t);
                case 1:
                {
                    dst[0] = (1.0 - y) * (1.0 - t);
                    dst[1] = -x * (1.0 - t);
                    return;
                }

                // f = x * y * (1.0 - t);
                case 2:
                {
                    dst[0] = y * (1.0 - t);
                    dst[1] = x * (1.0 - t);
                    return;
                }

                // f = (1.0 - x) * y * (1.0 - t);
                case 3:
                {
                    dst[0] = - y * (1.0 - t);
                    dst[1] = (1.0 - x) * (1.0 - t);
                    return;
                }

                // f = (1.0 - x) * (1.0 - y) * t;
                case 4:
                {
                    dst[0] = -(1.0 - y) * t;
                    dst[1] = -(1.0 - x) * t;
                    return;
                }

                // f = x * (1.0 - y) * t;
                case 5:
                {
                    dst[0] = (1.0 - y) * t;
                    dst[1] = -x * t;
                    return;
                }

                // f = x * y * t;
                case 6:
                {
                    dst[0] = y * t;
                    dst[1] = x * t;
                    return;
                }

                // f = (1.0 - x) * y * t;
                case 7:
                {
                    dst[0] = -y * t;
                    dst[1] = (1.0 - x) * t;
                    return;
                }

                default:
                {
                    dst[0] = 0.0;
                    dst[1] = 0.0;
                    return;
                }
            }

        }

        //space-time mixed derivative
        template<typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION static void grad_x_partial_t(const int i, const Point &p, Deriv &dst)
        {
            //project t coordinates to 0
            UTOPIA_DEVICE_ASSERT(dst.size() == 2);

            using Scalar = typename Traits<Point>::Scalar;

            const Scalar x = p[0];
            const Scalar y = p[1];
            const Scalar t = p[2];

            switch(i)
            {
                // f = (1.0 - x) * (1.0 - y) * (1.0 - t);
                // f_t = -(1.0 - x) * (1.0 - y)
                case 0:
                {
                    dst[0] = (1.0 - y);
                    dst[1] = (1.0 - x);
                    return;
                }

                // f = x * (1.0 - y) * (1.0 - t);
                // f_t = -x * (1.0 - y) 
                case 1:
                {
                    dst[0] = -(1.0 - y);
                    dst[1] = x;
                    return;
                }

                // f = x * y * (1.0 - t);
                // f_t = -x * y
                case 2:
                {
                    dst[0] = -y;
                    dst[1] = -x;
                    return;
                }

                // f = (1.0 - x) * y * (1.0 - t);
                // f_t = -(1.0 - x) * y 
                case 3:
                {
                    dst[0] = y;
                    dst[1] = -(1.0 - x);
                    return;
                }

                // f = (1.0 - x) * (1.0 - y) * t;
                // f_t = (1.0 - x) * (1.0 - y)
                case 4:
                {
                    dst[0] = -(1.0 - y);
                    dst[1] = -(1.0 - x);
                    return;
                }

                // f = x * (1.0 - y) * t;
                // f_t = x * (1.0 - y) 
                case 5:
                {
                    dst[0] = (1.0 - y);
                    dst[1] = -x;
                    return;
                }

                // f = x * y * t;
                // f_t = x * y
                case 6:
                {
                    dst[0] = y;
                    dst[1] = x;
                    return;
                }

                // f = (1.0 - x) * y * t;
                // f_t = (1.0 - x) * y 
                case 7:
                {
                    dst[0] = -y;
                    dst[1] = (1.0 - x);
                    return;
                }

                default:
                {
                    dst[0] = 0.0;
                    dst[1] = 0.0;
                    return;
                }
            }
        }

        template<typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION static void grad(const int i, const Point &p, Grad &g)
        {
            using Scalar = typename Traits<Point>::Scalar;

            const Scalar x = p[0];
            const Scalar y = p[1];
            const Scalar z = p[2];

            switch(i)
            {
                // f = (1.0 - x) * (1.0 - y) * (1.0 - z);
                case 0:
                {
                    g[0] = -(1.0 - y) * (1.0 - z);
                    g[1] = -(1.0 - x) * (1.0 - z);
                    g[2] = -(1.0 - x) * (1.0 - y);
                    return;
                }

                // f = x * (1.0 - y) * (1.0 - z);
                case 1:
                {
                    g[0] = (1.0 - y) * (1.0 - z);
                    g[1] = -x * (1.0 - z);
                    g[2] = -x * (1.0 - y);
                    return;
                }

                // f = x * y * (1.0 - z);
                case 2:
                {
                    g[0] = y * (1.0 - z);
                    g[1] = x * (1.0 - z);
                    g[2] = -x * y;
                    return;
                }

                // f = (1.0 - x) * y * (1.0 - z);
                case 3:
                {
                    g[0] = - y * (1.0 - z);
                    g[1] = (1.0 - x) * (1.0 - z);
                    g[2] = -(1.0 - x) * y;
                    return;
                }

                // f = (1.0 - x) * (1.0 - y) * z;
                case 4:
                {
                    g[0] = -(1.0 - y) * z;
                    g[1] = -(1.0 - x) * z;
                    g[2] = (1.0 - x) * (1.0 - y);
                    return;
                }

                // f = x * (1.0 - y) * z;
                case 5:
                {
                    g[0] = (1.0 - y) * z;
                    g[1] = -x * z;
                    g[2] = x * (1.0 - y);
                    return;
                }

                // f = x * y * z;
                case 6:
                {
                    g[0] = y * z;
                    g[1] = x * z;
                    g[2] = x * y;
                    return;
                }

                // f = (1.0 - x) * y * z;
                case 7:
                {
                    g[0] = -y * z;
                    g[1] = (1.0 - x) * z;
                    g[2] = (1.0 - x) * y;
                    return;
                }

                default:
                {
                    g[0] = 0.0;
                    g[1] = 0.0;
                    g[2] = 0.0;
                    return;
                }
            }
        }

    };

    template<typename Scalar_>
    class UniformHex8 : public Elem {
    public:
        using Scalar = Scalar_;
        using MemType = Uniform<>;
        static const int Dim = 3;
        static const int NNodes = 8;
        static const int NFunctions = NNodes;

        using Point     = utopia::StaticVector<Scalar, Dim>;
        using GradValue = utopia::StaticVector<Scalar, Dim>;
        using STGradX   = utopia::StaticVector<Scalar, Dim-1>;
        using FunValue  = Scalar;

        template<typename Point>
        UTOPIA_INLINE_FUNCTION static auto fun(const int i, const Point &p) -> decltype(RefHex8::fun(i, p))
        {
          return RefHex8::fun(i, p);
        }

        template<typename Point>
        UTOPIA_INLINE_FUNCTION void node(const std::size_t &i, Point &p) const
        {
            p[0] = translation_[0];
            p[1] = translation_[1];
            p[2] = translation_[2];

            switch(i) {
                case 0: { return; }
                case 1: { p[0] += h_[0]; return; }
                case 2: { p[0] += h_[0]; p[1] += h_[1]; return; }
                case 3: { p[1] += h_[1]; return; }
                case 4: { p[2] += h_[2]; return; }
                case 5: { p[2] += h_[2]; p[0] += h_[0]; return; }
                case 6: { p[2] += h_[2]; p[0] += h_[0]; p[1] += h_[1]; return; }
                case 7: { p[2] += h_[2]; p[1] += h_[1]; return; }
                default: {
                    UTOPIA_DEVICE_ASSERT(false);
                    return;
                }
            }
        }

        template<typename Point, typename Grad>
        UTOPIA_INLINE_FUNCTION void grad(const int i, const Point &p, Grad &g) const
        {
            RefHex8::grad(i, p, g);
            g[0] /= h_[0];
            g[1] /= h_[1];
            g[2] /= h_[2];
        }


        //space-time spatial gradient
        template<typename Point>
        UTOPIA_INLINE_FUNCTION auto partial_t(const int i, const Point &p) -> typename Traits<Point>::Scalar
        {
            return RefHex8::partial_t(i, p) / h_[2];
        }

        //space-time spatial gradient
        template<typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION void grad_x(const int i, const Point &p, Deriv &dst)
        {
            RefHex8::grad_x(i, p, dst);
            dst[0] /= h_[0];
            dst[1] /= h_[1];
        }

        ///space-time mixed derivative \nabla_x \partial_t \phi(x, t)
        template<typename Point, typename Deriv>
        UTOPIA_INLINE_FUNCTION void grad_x_partial_t(const int i, const Point &p, Deriv &dst)
        {
            RefHex8::grad_x_partial_t(i, p, dst);
            dst[0] /= (h_[0]*h_[2]);
            dst[1] /= (h_[1]*h_[2]);
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
            UTOPIA_DEVICE_ASSERT(h_[0]*h_[1]*h_[2] > 0.0);
            return h_[0]*h_[1]*h_[2];
        }

        template<typename RefPoint, typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void point(const RefPoint &in, PhysicalPoint &out) const
        {
            out[0] = in[0] * h_[0] + translation_[0];
            out[1] = in[1] * h_[1] + translation_[1];
            out[2] = in[2] * h_[2] + translation_[2];
        }

        template<typename PhysicalPoint>
        UTOPIA_INLINE_FUNCTION void centroid(PhysicalPoint &out) const
        {
            for(int i = 0; i < Dim; ++i) {
                out[i] = translation_[i] + h_[i]/2.0;
            }
        }

        UTOPIA_INLINE_FUNCTION UniformHex8()
        {
            h_[0] = 0.0;
            h_[1] = 0.0;
            h_[2] = 0.0;
        }

        template<class H>
        UTOPIA_INLINE_FUNCTION void set(
            const Point &translation,
            const H &h)
        {
            translation_(0) = translation(0);
            translation_(1) = translation(1);
            translation_(2) = translation(2);

            h_[0] = h[0];
            h_[1] = h[1];
            h_[2] = h[2];
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
        Scalar h_[3];
        Point translation_;
    };

}

#endif //UTOPIA_UNIFORM_HEX8_HPP
