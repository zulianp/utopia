#ifndef UTOPIA_SIMD_QUADRATURE_DB_HPP
#define UTOPIA_SIMD_QUADRATURE_DB_HPP

#include "utopia_simd_Quadrature_v2.hpp"

#include "utopia_MultiVariateElement.hpp"
// #include "utopia_Tet4.hpp"
#include "utopia_Tri3.hpp"
#include "utopia_UniformHex8.hpp"
#include "utopia_UniformQuad4.hpp"

namespace utopia {
    namespace simd_v2 {
        template <typename...>
        class QuadratureDB {};

        template <typename T, typename QT>
        class QuadratureDB<UniformHex8<T>, QT> {
        public:
            template <class Q>
            static bool get(const int order, Q &q) {
                return Gauss<QT>::Hex::get(order, q);
            }
        };

        template <typename T, typename QT>
        class QuadratureDB<UniformQuad4<T>, QT> {
        public:
            template <class Q>
            static bool get(const int order, Q &q) {
                return Gauss<QT>::Quad::get(order, q);
            }
        };

        template <typename T, typename QT>
        class QuadratureDB<Tri3<T>, QT> {
        public:
            template <class Q>
            static bool get(const int order, Q &q) {
                return Gauss<QT>::Tri::get(order, q);
            }
        };

        // template <typename T, typename QT>
        // class QuadratureDB<Tet4<T>, QT> {
        // public:
        //     template <class Q>
        //     static bool get(const int order, Q &q) {
        //         return Gauss<QT>::Tet::get(order, q);
        //     }
        // };

        template <typename Elem, int NVar, typename QT>
        class QuadratureDB<MultiVariateElem<Elem, NVar>, QT> : public QuadratureDB<Elem, QT> {};
    }  // namespace simd_v2
}  // namespace utopia

#endif  // UTOPIA_SIMD_QUADRATURE_DB_HPP
