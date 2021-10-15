

#include "utopia_Main.hpp"

#include "utopia_intrepid2.hpp"
#include "utopia_stk.hpp"
#include "utopia_stk_Commons.hpp"
#include "utopia_stk_intrepid2.hpp"

#include "utopia_kokkos_SamplePoints.hpp"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace utopia {

    template <class FunctionSpace>
    class SamplePointsApp {
    public:
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;
        using Field_t = utopia::Field<FunctionSpace>;
        using Intrepid2FE_t = utopia::intrepid2::FE<Scalar_t>;
        using DynRankView_t = typename Intrepid2FE_t::DynRankView;
        using Intrepid2Field_t = utopia::kokkos::Field<Intrepid2FE_t>;
        using SamplePoints_t = utopia::kokkos::SamplePoints<Intrepid2FE_t>;
        using HPolytope_t = utopia::kokkos::HPolytope<DynRankView_t>;

        static void create_normals_and_distances(const Mesh_t &mesh, DynRankView_t &normals, DynRankView_t &distances) {
            // TODO
            auto &bulk_data = mesh.bulk_data();
            auto &meta_data = mesh.meta_data();

            auto &&element_buckets = utopia::stk::local_elements(bulk_data);

            for (auto &&b_ptr : element_buckets) {
                auto &b = *b_ptr;

                const Size_t length = b.size();

                for (Size_t k = 0; k < length; ++k) {
                    auto entity = b[k];

                    bool has_face = false;

                    utopia::out() << "elem " << entity.local_offset() << ":\n";

                    for (auto f_ptr = bulk_data.begin_faces(entity); f_ptr != bulk_data.end_faces(entity); ++f_ptr) {
                        auto face_entity = *f_ptr;
                        utopia::out() << face_entity.local_offset() << " ";
                        has_face = true;
                    }

                    if (has_face) {
                        utopia::out() << "\n";
                    }
                }
            }
        }

        // void create_normals_and_distances(const Intrepid2FE_t &fe, DynRankView_t &normals, DynRankView_t &distances)
        // {
        //     Size_t n_cells = fe.n_cells();
        //     Size_t n_nodes_x_elem = fe.n_nodes_x_element();
        //     Size_t n_faces_x_element = fe.n_nodes_x_element();
        //     const int dim = fe.spatial_dimension();

        //     // Only simplex for now
        //     assert(dim + 1 == int(n_nodes_x_elem));

        //     auto points = fe.points();

        //     normals = DynRankView_t("normals", n_cells, n_faces_x_element, dim);
        //     distance = DynRankView_t("distance", n_cells, n_faces_x_element);

        //     ::Kokkos::parallel_for(
        //         fe.cell_range(), UTOPIA_LAMBDA(const Size_t cell) {
        //             Scalar_t u[3][3];

        //             // for(int s = 0; s < dim+1; ++s) {

        //             // }

        //             switch (dim) {
        //                 case 2: {
        //                     // Perp

        //                     break;
        //                 }
        //                 case 3: {
        //                     // Cross
        //                     break;
        //                 }
        //                 default: {
        //                     assert(false);
        //                     break;
        //                 }
        //             }
        //         });
        // }

        void run(Input &in) {
            FunctionSpace space;
            Field_t field;

            in.require("space", [&](Input &node) { space.read_with_state(node, field); });

            auto fe = std::make_shared<Intrepid2FE_t>();

            DynRankView_t normals, distances;
            create_normals_and_distances(space.mesh(), normals, distances);
            HPolytope_t h_poly(normals, distances);
            SamplePoints_t sample(fe, h_poly);
            sample.read(in);

            DynRankView_t values;
            Intrepid2Field_t fe_field(fe);
            convert_field(field, fe_field);

            sample.sample(fe_field, values);
        }
    };

}  // namespace utopia

void stk_sample_points(utopia::Input &in) {
    utopia::SamplePointsApp<utopia::stk::FunctionSpace> app;
    app.run(in);
}

UTOPIA_REGISTER_APP(stk_sample_points);
