#ifndef UTOPIA_KOKKOS_SAMPLE_POINTS_HPP
#define UTOPIA_KOKKOS_SAMPLE_POINTS_HPP

#include "utopia_Input.hpp"

namespace utopia {
    namespace kokkos {
        template <class FE_>
        class SamplePoints : public Configurable {
        public:
            using FE = FE_;
            using DynRankView = typename FE::DynRankView;
            using Scalar = typename FE::Scalar;

            using HPolytope = utopia::kokkos::HPolytope<DynRankView>;
            using Points = DynRankView;
            using HostPoints = typename Points::HostMirror;
            static constexpr int MAX_DIM = 4;

            void read(Input &in) override {
                int n_points = 0;
                in.get("points", [&](Input &array_node) { in.get_all([&](Input &) { ++n_points; }) });

                if (n_points == 0) {
                    assert(false);
                    return;
                }

                const int dim = fe_->spatial_dimension();
                points_ = Points("points", n_points, dim);
                HostPoints host_points = Kokkos::create_mirror_view(points_);

                in.get("points", [&](Input &array_node) {
                    int idx = 0;
                    in.get_all([&](Input &node) {
                        Scalar p[MAX_DIM] = {0., 0., 0., 0};

                        node.get("x", p[0]);
                        node.get("y", p[1]);
                        node.get("z", p[2]);
                        node.get("t", p[3]);

                        for (int d = 0; d < dim; ++d) {
                            host_points(idx, d) = p[d];
                        }
                    })
                });

                ::Kokkos::deep_copy(points_, host_points);
            }

            void sample(const Field<FE> field, DynRankView &values) const {
                assert(field.is_coefficient());
                // TODO
                values = DynRankView("SamplePoints_values", points_.extent(0), field.tensor_size());
            }

            void set_points(const Points &points) { points_ = points; }

            SamplePoints(const std::shared_ptr<FE> &fe, const HPolytope &polytope) : fe_(fe), polytope_(polytope) {}

        private:
            std::shared_ptr<FE> fe_;
            HPolytope polytope_;
            Points points_;
            Path output_path_{"out.csv"};
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_SAMPLE_POINTS_HPP