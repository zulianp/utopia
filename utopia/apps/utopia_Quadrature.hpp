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


    template<typename Scalar_, int Order, int Dim, int NPoints>
    class Hex8Quadrature {};

    template<typename Scalar_>
    class Hex8Quadrature<Scalar_, 2, 3, 27> {
    public:
        static const int Order   = 2;
        static const int Dim     = 3;
        static const int NPoints = 27;

        template<class Points, class Weights>
        UTOPIA_INLINE_FUNCTION static void get(Points &points, Weights &weights)
        {
            using PA = utopia::VectorAccessor<Points>;
            using WA = utopia::Accessor<Weights>;

            PA::set(points, 0, 0.112701665379258, 0.112701665379258, 0.112701665379258);
            PA::set(points, 1, 0.500000000000000, 0.112701665379258, 0.112701665379258);
            PA::set(points, 2, 0.887298334620742, 0.112701665379258, 0.112701665379258);
            PA::set(points, 3, 0.112701665379258, 0.500000000000000, 0.112701665379258);
            PA::set(points, 4, 0.500000000000000, 0.500000000000000, 0.112701665379258);
            PA::set(points, 5, 0.887298334620742, 0.500000000000000, 0.112701665379258);
            PA::set(points, 6, 0.112701665379258, 0.887298334620742, 0.112701665379258);
            PA::set(points, 7, 0.500000000000000, 0.887298334620742, 0.112701665379258);
            PA::set(points, 8, 0.887298334620742, 0.887298334620742, 0.112701665379258);
            PA::set(points, 9, 0.112701665379258, 0.112701665379258, 0.500000000000000);
            PA::set(points, 10, 0.500000000000000, 0.112701665379258, 0.500000000000000);
            PA::set(points, 11, 0.887298334620742, 0.112701665379258, 0.500000000000000);
            PA::set(points, 12, 0.112701665379258, 0.500000000000000, 0.500000000000000);
            PA::set(points, 13, 0.500000000000000, 0.500000000000000, 0.500000000000000);
            PA::set(points, 14, 0.887298334620742, 0.500000000000000, 0.500000000000000);
            PA::set(points, 15, 0.112701665379258, 0.887298334620742, 0.500000000000000);
            PA::set(points, 16, 0.500000000000000, 0.887298334620742, 0.500000000000000);
            PA::set(points, 17, 0.887298334620742, 0.887298334620742, 0.500000000000000);
            PA::set(points, 18, 0.112701665379258, 0.112701665379258, 0.887298334620742);
            PA::set(points, 19, 0.500000000000000, 0.112701665379258, 0.887298334620742);
            PA::set(points, 20, 0.887298334620742, 0.112701665379258, 0.887298334620742);
            PA::set(points, 21, 0.112701665379258, 0.500000000000000, 0.887298334620742);
            PA::set(points, 22, 0.500000000000000, 0.500000000000000, 0.887298334620742);
            PA::set(points, 23, 0.887298334620742, 0.500000000000000, 0.887298334620742);
            PA::set(points, 24, 0.112701665379258, 0.887298334620742, 0.887298334620742);
            PA::set(points, 25, 0.500000000000000, 0.887298334620742, 0.887298334620742);
            PA::set(points, 26, 0.887298334620742, 0.887298334620742, 0.88729833462074);


            WA::set(weights, 0, 0.021433470507545);
            WA::set(weights, 1, 0.034293552812071);
            WA::set(weights, 2, 0.021433470507545);
            WA::set(weights, 3, 0.034293552812071);
            WA::set(weights, 4, 0.054869684499314);
            WA::set(weights, 5, 0.034293552812071);
            WA::set(weights, 6, 0.021433470507545);
            WA::set(weights, 7, 0.034293552812071);
            WA::set(weights, 8, 0.021433470507545);
            WA::set(weights, 9, 0.034293552812071);
            WA::set(weights, 10, 0.054869684499314);
            WA::set(weights, 11, 0.034293552812071);
            WA::set(weights, 12, 0.054869684499314);
            WA::set(weights, 13, 0.087791495198903);
            WA::set(weights, 14, 0.054869684499314);
            WA::set(weights, 15, 0.034293552812071);
            WA::set(weights, 16, 0.054869684499314);
            WA::set(weights, 17, 0.034293552812071);
            WA::set(weights, 18, 0.021433470507545);
            WA::set(weights, 19, 0.034293552812071);
            WA::set(weights, 20, 0.021433470507545);
            WA::set(weights, 21, 0.034293552812071);
            WA::set(weights, 22, 0.054869684499314);
            WA::set(weights, 23, 0.034293552812071);
            WA::set(weights, 24, 0.021433470507545);
            WA::set(weights, 25, 0.034293552812071);
            WA::set(weights, 26, 0.021433470507545);
        }
    };


}

#endif //UTOPIA_QUADRATURE_HPP
