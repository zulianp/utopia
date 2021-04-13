#ifndef UTOPIA_SIMD_FE_TRAITS_HPP
#define UTOPIA_SIMD_FE_TRAITS_HPP

#include "utopia_simd_Traits.hpp"

#include "utopia_simd_Matrix.hpp"
#include "utopia_simd_Vector.hpp"

#include "utopia_simd_QuadratureDB.hpp"
#include "utopia_simd_Quadrature_v2.hpp"

namespace utopia {

    namespace simd_v2 {

        template <class Elem>
        struct FETraits {
            static const int Dim = Elem::Dim;
            using T = typename Elem::Scalar;
            using Point = simd_v2::Vector<T, Dim>;
            using Scalar = T;
            using SIMDType = typename Traits<Point>::Scalar;

            using FunValue = SIMDType;
            using GradValue = Point;
            // using STGradX = simd_v2::Vector<T, Dim - 1>;
        };

        template <class Elem, int NVar>
        struct FETraits<MultiVariateElem<Elem, NVar>> {
            static const int Dim = Elem::Dim;
            using T = typename Elem::Scalar;
            using Point = simd_v2::Vector<T, Dim>;
            using Scalar = typename Traits<T>::Scalar;
            using FunValue = simd_v2::Vector<T, NVar>;
            using GradValue = simd_v2::Matrix<T, NVar, Dim>;
            // using STGradX = simd_v2::Vector<T, Dim - 1>;
        };

        template <class Elem>
        struct FETraits<MultiVariateElem<Elem, 1>> : FETraits<Elem> {};

    }  // namespace simd_v2

}  // namespace utopia

#endif  // UTOPIA_SIMD_FE_TRAITS_HPP
