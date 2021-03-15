#include "utopia_Main.hpp"
#include "utopia_fe_config.hpp"

#ifdef UTOPIA_WITH_LIBMESH
#include "MainOld.hpp"
#include "utopia_libmesh_Library.hpp"
#endif

int main(const int argc, char *argv[]) {
#ifdef UTOPIA_WITH_LIBMESH
    // utopia::Utopia::instance().add_library(utopia::make_unique<utopia::LibMeshLibrary>());
    return MainOld(argc, argv);
#endif  // UTOPIA_WITH_LIBMESH

    return UTOPIA_MAIN(argc, argv);
}
