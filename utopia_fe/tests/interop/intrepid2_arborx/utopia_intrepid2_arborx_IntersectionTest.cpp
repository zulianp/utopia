#include "utopia_Testing.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_ui.hpp"

#include "utopia_intrepid2_arborx_Intersector.hpp"

using namespace utopia;

using FE = utopia::intrepid2::FE<double>;

static std::shared_ptr<FE> make_ref_tet() {
    double tet[1 * 4 * 3] = {
        // p0
        0.0,
        0.0,
        0.0,
        // p1
        1.0,
        0.0,
        0.0,
        // p2
        0.0,
        1.0,
        0.0,
        // p3
        0.0,
        0.0,
        1.0,
    };

    int n_cells = 1;
    int n_nodes = 4;
    int n_dims = 3;

    FE::DynRankView device_cell_points("cell_points", n_cells, n_nodes, n_dims);
    FE::DynRankView::HostMirror cell_points = ::Kokkos::create_mirror_view(device_cell_points);

    for (int c = 0; c < n_cells; ++c) {
        for (int n = 0; n < n_nodes; ++n) {
            for (int d = 0; d < n_dims; ++d) {
                cell_points(c, n, d) = tet[c * (n_nodes * n_dims) + n * n_dims + d];
            }
        }
    }

    ::Kokkos::deep_copy(device_cell_points, cell_points);

    shards::CellTopology cell_type = shards::getCellTopologyData<shards::Tetrahedron<>>();

    auto fe_ptr = std::make_shared<FE>();
    fe_ptr->init(cell_type, device_cell_points);
    return fe_ptr;
}

void intersect_simple() {
    auto fe_ptr = make_ref_tet();

    Intrepid2ArborXIntersector<FE::DynRankView> isect;

    isect.detect(fe_ptr->cell_nodes, fe_ptr->cell_nodes);
}

void intrepid2arborx() { UTOPIA_RUN_TEST(intersect_simple); }

UTOPIA_REGISTER_TEST_FUNCTION(intrepid2arborx);
