#include "utopia.hpp"
#include "utopia_Testing.hpp"

// FIXME remove mre
#include <mpi.h>
#include "libmesh/libmesh.h"

int main(int argc, char *argv[]) {
    using namespace utopia;

    // For debugginh with ddt
    MPI_Init(&argc, &argv);
    PETSC_COMM_WORLD = MPI_COMM_WORLD;

    Utopia::Init(argc, argv);

    {
        // FIXME remove me
        libMesh::LibMeshInit init(argc, argv, PETSC_COMM_WORLD);

        Utopia::Init(argc, argv);
        UTOPIA_TRACE_REGION_BEGIN("main");

        {
            TestRunner app_runner;
            app_runner.run(argc, argv);
        }
    }

    return Utopia::Finalize();
}
