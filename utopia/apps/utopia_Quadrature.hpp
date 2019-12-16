#ifndef UTOPIA_QUADRATURE_HPP
#define UTOPIA_QUADRATURE_HPP

#include "utopia_Accessor.hpp"

namespace utopia {

    template<typename T, int Order, int Dim = T::Dim, typename ...>
    class Quadrature {};

    template<typename Scalar_, int Order, int Dim, int NPoints>
    class Quad4Quadrature {};

    template<typename Scalar_>
    class Quad4Quadrature<Scalar_, 2, 2, 6> {
    public:
        static const int Order   = 2;
        static const int Dim     = 2;
        static const int NPoints = 6;

        template<class Points, class Weights>
        UTOPIA_INLINE_FUNCTION static void get(Points &points, Weights &weights)
        {
            using PA = utopia::Accessor<Points>;
            using WA = utopia::Accessor<Weights>;

            PA::set(points, 0, 0, 0.5);
            PA::set(points, 0, 1, 0.5);
            PA::set(points, 1, 0, 0.9304589153964795245728880523899),
            PA::set(points, 1, 1, 0.5);
            PA::set(points, 2, 0, 0.72780186391809642112479237299488);
            PA::set(points, 2, 1, 0.074042673347699754349082179816666);
            PA::set(points, 3, 0, 0.72780186391809642112479237299488);
            PA::set(points, 3, 1, 0.92595732665230024565091782018333);
            PA::set(points, 4, 0, 0.13418502421343273531598225407969);
            PA::set(points, 4, 1, 0.18454360551162298687829339850317);
            PA::set(points, 5, 0, 0.13418502421343273531598225407969);
            PA::set(points, 5, 1, 0.81545639448837701312170660149683);

            WA::set(weights, 0, 0.28571428571428571428571428571428);
            WA::set(weights, 1, 0.10989010989010989010989010989011);
            WA::set(weights, 2, 0.14151805175188302631601261486295);
            WA::set(weights, 3, 0.14151805175188302631601261486295);
            WA::set(weights, 4, 0.16067975044591917148618518733485);
            WA::set(weights, 5, 0.16067975044591917148618518733485);
        }
    };
}

#endif //UTOPIA_QUADRATURE_HPP
