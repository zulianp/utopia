#include "utopia_Testing.hpp"

#include "utopia_intrepid2_LaplaceOperator.hpp"
#include "utopia_ui.hpp"

using namespace utopia;

using HGradFE = utopia::intrepid2::HGradFE<double>;

static std::shared_ptr<HGradFE> make_ref_tet() {
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

    HGradFE::DynRankView cell_points("cell_points", n_cells, n_nodes, n_dims);

    for (int c = 0; c < n_cells; ++c) {
        for (int n = 0; n < n_nodes; ++n) {
            for (int d = 0; d < n_dims; ++d) {
                cell_points(c, n, d) = tet[c * (n_nodes * n_dims) + n * n_dims + d];
            }
        }
    }

    shards::CellTopology cell_type = shards::getCellTopologyData<shards::Tetrahedron<>>();

    auto fe_ptr = std::make_shared<HGradFE>();
    fe_ptr->init(cell_type, cell_points);
    return fe_ptr;
}

void intrepid2_basis_functions() {
    auto fe_ptr = make_ref_tet();

    auto actual = fe_ptr->measure(0, 0);
    utopia_test_assert(approxeq(1.0 / 6.0, actual, 1e-10));
}

void intrepid2_laplace_operator() {
    LaplaceOperator<double> lapl{1.0};
    auto fe_ptr = make_ref_tet();

    intrepid2::Assemble<LaplaceOperator<double>> assembler(lapl, fe_ptr);
    assembler.init();
    assembler.describe(std::cout);
}

void uintrepid2() {
    UTOPIA_RUN_TEST(intrepid2_basis_functions);
    UTOPIA_RUN_TEST(intrepid2_laplace_operator);
}

UTOPIA_REGISTER_TEST_FUNCTION(uintrepid2);
