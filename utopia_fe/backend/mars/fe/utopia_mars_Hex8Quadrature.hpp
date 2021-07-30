#ifndef UTOPIA_MARS_HEX_8_QUADRATURE_HPP
#define UTOPIA_MARS_HEX_8_QUADRATURE_HPP

#include "utopia_Accessor.hpp"

namespace utopia {

    template <typename Scalar, int Order, int Dim, int NPoints>
    class Hex8Quadrature {};

    template <typename Scalar>
    class Hex8Quadrature<Scalar, 2, 3, 27> {
    public:
        static const int Order = 2;
        static const int Dim = 3;
        static const int NPoints = 27;

        template <class Array>
        UTOPIA_INLINE_FUNCTION static void set(Array &a, const int i, const Scalar &value) {
            a(i) = value;
        }

        template <class Array>
        UTOPIA_INLINE_FUNCTION static void set(Array &a,
                                               const int i,
                                               const Scalar &x,
                                               const Scalar &y,
                                               const Scalar &z) {
            a(i, 0) = x;
            a(i, 1) = y;
            a(i, 2) = z;
        }

        template <class Points, class Weights>
        UTOPIA_INLINE_FUNCTION static void get(Points &points, Weights &weights) {
            set(points, 0, 0.112701665379258, 0.112701665379258, 0.112701665379258);
            set(points, 1, 0.500000000000000, 0.112701665379258, 0.112701665379258);
            set(points, 2, 0.887298334620742, 0.112701665379258, 0.112701665379258);
            set(points, 3, 0.112701665379258, 0.500000000000000, 0.112701665379258);
            set(points, 4, 0.500000000000000, 0.500000000000000, 0.112701665379258);
            set(points, 5, 0.887298334620742, 0.500000000000000, 0.112701665379258);
            set(points, 6, 0.112701665379258, 0.887298334620742, 0.112701665379258);
            set(points, 7, 0.500000000000000, 0.887298334620742, 0.112701665379258);
            set(points, 8, 0.887298334620742, 0.887298334620742, 0.112701665379258);
            set(points, 9, 0.112701665379258, 0.112701665379258, 0.500000000000000);
            set(points, 10, 0.500000000000000, 0.112701665379258, 0.500000000000000);
            set(points, 11, 0.887298334620742, 0.112701665379258, 0.500000000000000);
            set(points, 12, 0.112701665379258, 0.500000000000000, 0.500000000000000);
            set(points, 13, 0.500000000000000, 0.500000000000000, 0.500000000000000);
            set(points, 14, 0.887298334620742, 0.500000000000000, 0.500000000000000);
            set(points, 15, 0.112701665379258, 0.887298334620742, 0.500000000000000);
            set(points, 16, 0.500000000000000, 0.887298334620742, 0.500000000000000);
            set(points, 17, 0.887298334620742, 0.887298334620742, 0.500000000000000);
            set(points, 18, 0.112701665379258, 0.112701665379258, 0.887298334620742);
            set(points, 19, 0.500000000000000, 0.112701665379258, 0.887298334620742);
            set(points, 20, 0.887298334620742, 0.112701665379258, 0.887298334620742);
            set(points, 21, 0.112701665379258, 0.500000000000000, 0.887298334620742);
            set(points, 22, 0.500000000000000, 0.500000000000000, 0.887298334620742);
            set(points, 23, 0.887298334620742, 0.500000000000000, 0.887298334620742);
            set(points, 24, 0.112701665379258, 0.887298334620742, 0.887298334620742);
            set(points, 25, 0.500000000000000, 0.887298334620742, 0.887298334620742);
            set(points, 26, 0.887298334620742, 0.887298334620742, 0.88729833462074);

            set(weights, 0, 0.021433470507545);
            set(weights, 1, 0.034293552812071);
            set(weights, 2, 0.021433470507545);
            set(weights, 3, 0.034293552812071);
            set(weights, 4, 0.054869684499314);
            set(weights, 5, 0.034293552812071);
            set(weights, 6, 0.021433470507545);
            set(weights, 7, 0.034293552812071);
            set(weights, 8, 0.021433470507545);
            set(weights, 9, 0.034293552812071);
            set(weights, 10, 0.054869684499314);
            set(weights, 11, 0.034293552812071);
            set(weights, 12, 0.054869684499314);
            set(weights, 13, 0.087791495198903);
            set(weights, 14, 0.054869684499314);
            set(weights, 15, 0.034293552812071);
            set(weights, 16, 0.054869684499314);
            set(weights, 17, 0.034293552812071);
            set(weights, 18, 0.021433470507545);
            set(weights, 19, 0.034293552812071);
            set(weights, 20, 0.021433470507545);
            set(weights, 21, 0.034293552812071);
            set(weights, 22, 0.054869684499314);
            set(weights, 23, 0.034293552812071);
            set(weights, 24, 0.021433470507545);
            set(weights, 25, 0.034293552812071);
            set(weights, 26, 0.021433470507545);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_MARS_HEX_8_QUADRATURE_HPP
