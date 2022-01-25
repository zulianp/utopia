#include "utopia_Testing.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_kokkos_LaplaceOperator.hpp"
#include "utopia_ui.hpp"

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

void intrepid2_basis_functions() {
    auto fe_ptr = make_ref_tet();

    FE::DynRankView::HostMirror host_measure = ::Kokkos::create_mirror_view(fe_ptr->measure());
    ::Kokkos::deep_copy(host_measure, fe_ptr->measure());

    auto actual = host_measure(0, 0);

    utopia_test_assert(approxeq(1.0 / 6.0, actual, 1e-10));
}

void intrepid2_laplace_operator() {
    kokkos::LaplaceOperator<FE>::Params lapl{1.0};
    auto fe_ptr = make_ref_tet();

    kokkos::LaplaceOperator<FE> assembler(fe_ptr, lapl);
    assembler.assemble_matrix();
    // assembler.describe(std::cout);

    // if (std::is_same<FE::ExecutionSpace, Kokkos::Cuda>::value) {
    //     utopia::out() << "Cuda\n";
    // } else if (std::is_same<FE::ExecutionSpace, Kokkos::OpenMP>::value) {
    //     utopia::out() << "OpenMP\n";
    // } else {
    //     utopia::out() << "What?\n";
    // }
}

void uintrepid2() {
    UTOPIA_RUN_TEST(intrepid2_basis_functions);
    UTOPIA_RUN_TEST(intrepid2_laplace_operator);
}

UTOPIA_REGISTER_TEST_FUNCTION(uintrepid2);
