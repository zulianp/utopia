#include "utopia.hpp"
#include "utopia_AppRunner.hpp"
#include "utopia_fe_Instance.hpp"
// FIXME remove mre
#include <mpi.h>
#include "libmesh/libmesh.h"

int main(int argc, char *argv[]) {
    using namespace utopia;

    UtopiaFE::Init(argc, argv);
    UTOPIA_TRACE_REGION_BEGIN("main");

    {
        AppRunner app_runner;
        app_runner.run(argc, argv);
    }

    UTOPIA_TRACE_REGION_END("main");
    return UtopiaFE::Finalize();
}
