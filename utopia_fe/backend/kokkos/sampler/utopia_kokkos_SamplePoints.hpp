#ifndef UTOPIA_KOKKOS_SAMPLE_POINTS_HPP
#define UTOPIA_KOKKOS_SAMPLE_POINTS_HPP

#include "utopia_Input.hpp"
#include "utopia_kokkos_HPolytope.hpp"
#include "utopia_kokkos_LinearSimplex.hpp"
#include "utopia_kokkos_PointHPolytopeIntersection.hpp"

namespace utopia {
    namespace kokkos {
        template <class FE_>
        class SamplePoints : public Configurable {
        public:
            using FE = FE_;
            using DynRankView = typename FE::DynRankView;
            using Scalar = typename FE::Scalar;
            using ExecutionSpace = typename FE::ExecutionSpace;

            using HPolytope = utopia::kokkos::HPolytope<DynRankView>;
            using Points = DynRankView;
#if KOKKOS_VERSION >= 50002
            using HostPoints = typename Points::host_mirror_type;
#else
            using HostPoints = typename Points::HostMirror;
#endif
            using Intersector = utopia::kokkos::PointHPolytopeIntersection<Points, HPolytope>;
            using Index = typename Intersector::IndexView;

            using Range1 = ::Kokkos::RangePolicy<ExecutionSpace>;
            using Range2 = ::Kokkos::MDRangePolicy<::Kokkos::Rank<2>, ExecutionSpace>;

            static constexpr int MAX_DIM = 4;

            void read(Input &in) override {
                int n_points = 0;
                in.get("points", [&](Input &array_node) { array_node.get_all([&](Input &) { ++n_points; }); });

                if (n_points == 0) {
                    assert(false);
                    return;
                }

                const int dim = fe_->spatial_dimension();
                points_ = Points("points", n_points, dim);
                HostPoints host_points = Kokkos::create_mirror_view(points_);

                in.get("points", [&](Input &array_node) {
                    int idx = 0;
                    array_node.get_all([&](Input &node) {
                        Scalar p[MAX_DIM] = {0., 0., 0., 0};

                        node.get("x", p[0]);
                        node.get("y", p[1]);
                        node.get("z", p[2]);
                        node.get("t", p[3]);

                        for (int d = 0; d < dim; ++d) {
                            host_points(idx, d) = p[d];
                        }
                    });
                });

                ::Kokkos::deep_copy(points_, host_points);
            }

            void sample(const Field<FE> &field, DynRankView &values) const {
                assert(field.is_coefficient());

                switch (fe_->spatial_dimension()) {
                    case 1: {
                        simplex_sample(Simplex<DynRankView, 1>(fe_->points()), field, values);
                        break;
                    }
                    case 2: {
                        simplex_sample(Simplex<DynRankView, 2>(fe_->points()), field, values);
                        break;
                    }

                    case 3: {
                        simplex_sample(Simplex<DynRankView, 3>(fe_->points()), field, values);
                        break;
                    }

                    case 4: {
                        simplex_sample(Simplex<DynRankView, 4>(fe_->points()), field, values);
                        break;
                    }
                    default: {
                        assert(false);
                        Utopia::Abort("Invalid dimension!");
                    }
                }
            }

            template <int Dim>
            void simplex_sample(const Simplex<DynRankView, Dim> &simplex,
                                const Field<FE> &field,
                                DynRankView &values) const {
                const SizeType n_points = points_.extent(0);
                const int tensor_size = field.tensor_size();

                values = DynRankView("SamplePoints_values", n_points, tensor_size);

                auto data = field.data();

                ::Kokkos::parallel_for(
                    Range1(0, points_.extent(0)), UTOPIA_LAMBDA(const SizeType i) {
                        SizeType idx = index_(i);
                        if (idx == -1) {
                            for (int t = 0; t < tensor_size; ++t) {
                                values(i, t) = 0.;
                            }
                        }

                        StaticVector<Scalar, Dim> p, p_ref;
                        for (int d = 0; d < Dim; ++d) {
                            p[d] = points_(i, d);
                        }

                        simplex.inverse_transform(idx, p, p_ref);

                        for (int t = 0; t < tensor_size; ++t) {
                            Scalar interpolated = 0.0;
                            for (int f = 0; f < Dim + 1; ++f) {
                                Scalar v = data(idx, f * tensor_size + t);
                                interpolated += v * simplex.fun(f, p_ref);
                            }

                            values(i, t) = interpolated;
                        }
                    });
            }

            void set_points(const Points &points) {
                points_ = points;
                find_point_locations();
            }

            SamplePoints(const std::shared_ptr<FE> &fe, const HPolytope &polytope) : fe_(fe), polytope_(polytope) {}

        private:
            std::shared_ptr<FE> fe_;
            HPolytope polytope_;
            Points points_;
            Index index_;
            // Path output_path_{"out.csv"};

            void find_point_locations() {
                Intersector intersector;
                intersector.intersect(points_, polytope_, index_);
            }
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_SAMPLE_POINTS_HPP