#include "utopia_Main.hpp"
#include "utopia_fe_config.hpp"

#ifdef UTOPIA_ENABLE_LIBMESH
#include "utopia_libmesh_Library.hpp"
#endif

#ifdef UTOPIA_ENABLE_MARS
#include "utopia_mars_Library.hpp"
#endif

int main(const int argc, char *argv[]) {
#ifdef UTOPIA_ENABLE_LIBMESH
    utopia::Utopia::instance().add_library(utopia::make_unique<utopia::LibMeshLibrary>());
#endif  // #ifdef UTOPIA_ENABLE_LIBMESH

    // #ifdef UTOPIA_ENABLE_MARS
    //     utopia::Utopia::instance().add_library(utopia::make_unique<utopia::MarsLibrary>());
    // #endif  // UTOPIA_ENABLE_MARS

    return UTOPIA_MAIN(argc, argv);
}
