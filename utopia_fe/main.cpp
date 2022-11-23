#include "utopia_Main.hpp"
#include "utopia_fe_config.hpp"

#ifdef UTOPIA_WITH_LIBMESH
#include "utopia_libmesh_Library.hpp"
#endif

#ifdef UTOPIA_WITH_MARS
#include "utopia_mars_Library.hpp"
#endif

int main(const int argc, char *argv[]) {
#ifdef UTOPIA_WITH_LIBMESH
    utopia::Utopia::instance().add_library(utopia::make_unique<utopia::LibMeshLibrary>());
#endif  //#ifdef UTOPIA_WITH_LIBMESH

    // #ifdef UTOPIA_WITH_MARS
    //     utopia::Utopia::instance().add_library(utopia::make_unique<utopia::MarsLibrary>());
    // #endif  // UTOPIA_WITH_MARS

    return UTOPIA_MAIN(argc, argv);
}
