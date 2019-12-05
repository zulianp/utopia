#ifndef UTOPIA_TRILINOS_STRUCTURED_GRID_HPP
#define UTOPIA_TRILINOS_STRUCTURED_GRID_HPP

#include "utopia_trilinos.hpp"
#include "utopia_trilinos_FECrsGraph.hpp"
#include "utopia_Views.hpp"

namespace utopia {

    class FE {};
    class FD {};
    class FV {};

    template<class Elem_, class MemType = typename Elem_::MemType, typename...>
    class Mesh {
    public:
        using Elem = Elem_;
    };

    template<class...>
    class FunctionSpace {};

    template<typename T, int Order, int Dim, int NPoints, typename ...>
    class Quadrature {};

    template<typename...>
    class Uniform {};

    template<typename...>
    class Varying {};

    class RefQuad4 {
    public:
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

    template<typename Scalar>
    class UniformQuad4 final {
    public:
        using MemType = Uniform<>;
        using DiscretizationType = FE;

        template<typename Point, typename Values>
        UTOPIA_INLINE_FUNCTION static void fun(const Point &p, Values &values)
        {
            RefQuad4::fun(p, values);
        }

        template<typename Point, typename Values>
        UTOPIA_INLINE_FUNCTION void grad(const Point &p, Values &values)
        {
            RefQuad4::grad(p, values);
            for(int i = 0; i < 4; ++i)
            {
                values(i, 0) /= h_[0];
                values(i, 1) /= h_[1];
            }
        }

        UTOPIA_INLINE_FUNCTION UniformQuad4(const Scalar &hx, const Scalar &hy)
        {
            h_[0] = hx;
            h_[1] = hy;
        }

    private:
        Scalar h_[2];
    };

    template<typename Scalar_>
    class Quadrature<UniformQuad4<Scalar_>, 2, 2, 6> {
    public:
        using Scalar = Scalar_;
        using PointView   = Kokkos::DualView<Scalar **>;
        using WeightView  = Kokkos::DualView<Scalar *>;


        static const int Order   = 2;
        static const int Dim     = 2;
        static const int NPoints = 6;

        void init()
        {
            auto host_points  = points_.view_host();
            auto host_weights = weights_.view_host();

            host_points(0, 0) = 0.5;
            host_points(0, 1) = 0.5;
            host_points(1, 0) = 0.9304589153964795245728880523899,
            host_points(1, 1) = 0.5;
            host_points(2, 0) = 0.72780186391809642112479237299488;
            host_points(2, 1) = 0.074042673347699754349082179816666;
            host_points(3, 0) = 0.72780186391809642112479237299488;
            host_points(3, 1) = 0.92595732665230024565091782018333;
            host_points(4, 0) = 0.13418502421343273531598225407969;
            host_points(4, 1) = 0.18454360551162298687829339850317;
            host_points(5, 0) = 0.13418502421343273531598225407969;
            host_points(5, 1) = 0.81545639448837701312170660149683;

            host_weights(0) = 0.28571428571428571428571428571428;
            host_weights(1) = 0.10989010989010989010989010989011;
            host_weights(2) = 0.14151805175188302631601261486295;
            host_weights(3) = 0.14151805175188302631601261486295;
            host_weights(4) = 0.16067975044591917148618518733485;
            host_weights(5) = 0.16067975044591917148618518733485;
        }

        Quadrature() :
            points_("Quadrature::points", NPoints, Dim),
            weights_("Quadrature::weights", NPoints)
        {}

    private:
        PointView points_;
        WeightView weights_;
    };

    template<class Mesh, class MemType = typename Mesh::Elem::MemType>
    class FEValues {};

    template<class Elem_, typename...Args>
    class Mesh< Elem_, Uniform<Args...>> {
    public:
        using Elem = Elem_;
    };

    template<class Elem, typename...Args>
    class FEValues<Mesh<Elem>, Uniform<Args...>> {
    public:
        using Scalar   = typename Elem::Scalar;
        using Mesh     = utopia::Mesh<Elem>;
        using FunView  = Kokkos::DualView<Scalar **>;
        using GradView = Kokkos::DualView<Scalar ***>;

        template<typename...GridArgs>
        void init(const FunctionSpace<Mesh> &space)
        {

        }

        // void
    };

    template<typename Scalar>
    using UniformQuad4Mesh = utopia::Mesh<UniformQuad4<Scalar>>;
}

#endif