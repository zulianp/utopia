
#include "utopia_Testing.hpp"

#include "utopia_MeshTest.hpp"

#include "utopia_petsc_Mesh.hpp"
#include "utopia_ui.hpp"

#include "utopia_RunParallelTest.hpp"

using namespace utopia;
using namespace utopia::petsc;

void unit_cube() {
    PetscCommunicator comm;
    const SizeType n = comm.size() * 2;

    StructuredGrid mesh(comm);
    mesh.unit_cube(n, n, n);

    mesh.describe(utopia::out().stream());

    const SizeType n_nodes = mesh.n_nodes();
    UTOPIA_TEST_EQ(n_nodes, ((n + 1) * (n + 1) * (n + 1)));
}

void upetsc() {
    // utopia::run_parallel_test<MeshTest<utopia::petsc::Mesh>>();
    unit_cube();
}

UTOPIA_REGISTER_TEST_FUNCTION(upetsc);
