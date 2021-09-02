#ifndef UTOPIA_KOKKOS_H_POLYTOPE_HPP
#define UTOPIA_KOKKOS_H_POLYTOPE_HPP

#include "utopia_Base.hpp"

#include "utopia_kokkos_Commons.hpp"

namespace utopia {
    namespace kokkos {
        template <class View>
        class HPolytope {
        public:
            using SizeType = typename Traits<View>::SizeType;
            using Scalar = typename Traits<View>::Scalar;

            UTOPIA_INLINE_FUNCTION Scalar distance(const SizeType idx, const SizeType side, const Scalar *point) const {
                Scalar ret = 0;
                for (int d = 0; d < dim; ++d) {
                    ret += point[d] * normals(idx, side, d);
                }

                ret -= distances(idx, side);
            }

            UTOPIA_INLINE_FUNCTION bool inside_half_space(const SizeType idx,
                                                          const SizeType side,
                                                          const Scalar *point) const {
                return distance(idx, side, point) <= tol;
            }

            UTOPIA_INLINE_FUNCTION bool contains(const SizeType idx, const Scalar *point) const {
                const SizeType n_sides = distances.extent(1);

                for (SizeType i = 0; i < n_sides; ++i) {
                    if (!inside_half_space(idx, i, point)) return false;
                }

                return true;
            }

            UTOPIA_INLINE_FUNCTION HPolytope(const View &normals, const View &distances, const Scalar &tol = 0.)
                : normals(normals), distances(distances), tol(tol), dim(distances.extent(2)) {}

            UTOPIA_INLINE_FUNCTION SizeType size() const { return normals.extent(0); }

            View normals;
            View distances;
            Scalar tol{0};
            int dim;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_H_POLYTOPE_HPP
