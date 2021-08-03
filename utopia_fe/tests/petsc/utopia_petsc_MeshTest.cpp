
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

void clone() {
    PetscCommunicator comm;
    const SizeType n = comm.size() * 2;

    StructuredGrid mesh(comm);
    mesh.unit_cube(n, n, n);

    auto c = mesh.clone();
    c->describe(utopia::out().stream());
}

void uniform_refine() {
    PetscCommunicator comm;
    const SizeType n = comm.size() * 2;

    StructuredGrid mesh(comm);
    mesh.unit_cube(n, n, n);

    auto c = mesh.uniform_refine();
    c->describe(utopia::out().stream());
}

void upetsc() {
    // utopia::run_parallel_test<MeshTest<utopia::petsc::Mesh>>();
    UTOPIA_RUN_TEST(unit_cube);
    UTOPIA_RUN_TEST(clone);
    UTOPIA_RUN_TEST(uniform_refine);
}

UTOPIA_REGISTER_TEST_FUNCTION(upetsc);
