#ifndef UTOPIA_KOKKOS_POINT_H_POLYTOPE_INTERSECTION_HPP
#define UTOPIA_KOKKOS_POINT_H_POLYTOPE_INTERSECTION_HPP

#include "utopia_kokkos_FEBase.hpp"

namespace utopia {
    namespace kokkos {
        template <class PointView, class HPolytopeView, class ExecutionSpace = typename PointView::execution_space>
        class PointHPolytopeIntersection {
        public:
            using SizeType = typename Traits<PointView>::SizeType;
            using Scalar = typename Traits<PointView>::Scalar;
            using IndexView = ::Kokkos::View<SizeType *, ExecutionSpace>;
            static constexpr int MAX_DIM = 4;

            /// Complexity is n*m
            void intersect(const PointView &points, const HPolytopeView &poly, IndexView &result) {
                using Range1 = ::Kokkos::RangePolicy<ExecutionSpace>;
                using Range2 = ::Kokkos::MDRangePolicy<::Kokkos::Rank<2>, ExecutionSpace>;

                const SizeType n_points = points.extent(0);
                const SizeType n_polytopes = poly.size();
                const int dim = points.extent(1);

                assert(dim < MAX_DIM);

                result = IndexView("PointHPolytopeIntersection_result", n_points);

                ::Kokkos::parallel_for(
                    Range1(0, n_points), UTOPIA_LAMBDA(const SizeType &idx) { result(idx) = -1; });

                ::Kokkos::parallel_for(
                    Range2({0, n_polytopes}, {0, n_points}),
                    UTOPIA_LAMBDA(const SizeType &idx_poly, const SizeType &idx_point) {
                        Scalar p[MAX_DIM];

                        for (int d = 0; d < dim; ++d) {
                            p[d] = points(idx_point, d);
                        }

                        if (poly.contains(idx_poly, p)) {
                            // Is this safe?
                            volatile SizeType *const r = &result(idx_point);
                            *r = idx_poly;
                        }
                    });
            }
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_POINT_H_POLYTOPE_INTERSECTION_HPP
