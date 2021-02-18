#ifndef UTOPIA_SIMD_QUADRATURE_DB_HPP
#define UTOPIA_SIMD_QUADRATURE_DB_HPP

#include "utopia_simd_Quadrature_v2.hpp"

#include "utopia_MultiVariateElement.hpp"
// #include "utopia_Tet4.hpp"
#include "utopia_Tri3.hpp"
#include "utopia_UniformHex8.hpp"
#include "utopia_UniformQuad4.hpp"

#include "utopia_simd_UniformHex8.hpp"
#include "utopia_simd_UniformQuad4.hpp"

namespace utopia {
    namespace simd_v2 {
        template <typename...>
        class QuadratureDB {};

        template <typename T>
        class QuadratureDB<utopia::UniformHex8<T>> {
        public:
            template <class Q>
            static bool get(const int order, Q &q) {
                return Gauss<T>::Hex::get(order, q);
            }
        };

        template <typename T>
        class QuadratureDB<utopia::UniformQuad4<T>> {
        public:
            template <class Q>
            static bool get(const int order, Q &q) {
                return Gauss<T>::Quad::get(order, q);
            }
        };

        // template <typename T>
        // class QuadratureDB<simd_v2::UniformHex8<T>> {
        // public:
        //     template <class Q>
        //     static bool get(const int order, Q &q) {
        //         return Gauss<T>::Hex::get(order, q);
        //     }
        // };

        // template <typename T>
        // class QuadratureDB<simd_v2::UniformQuad4<T>> {
        // public:
        //     template <class Q>
        //     static bool get(const int order, Q &q) {
        //         return Gauss<T>::Quad::get(order, q);
        //     }
        // };

        template <typename T>
        class QuadratureDB<Tri3<T>> {
        public:
            template <class Q>
            static bool get(const int order, Q &q) {
                return Gauss<T>::Tri::get(order, q);
            }
        };

        // template <typename T>
        // class QuadratureDB<Tet4<T>> {
        // public:
        //     template <class Q>
        //     static bool get(const int order, Q &q) {
        //         return Gauss<T>::Tet::get(order, q);
        //     }
        // };

        template <typename Elem, int NVar>
        class QuadratureDB<MultiVariateElem<Elem, NVar>> : public QuadratureDB<Elem> {};
    }  // namespace simd_v2
}  // namespace utopia

#endif  // UTOPIA_SIMD_QUADRATURE_DB_HPP
