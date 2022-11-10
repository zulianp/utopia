#include "utopia_Main.hpp"
#include "utopia_fe_config.hpp"

#ifdef UTOPIA_WITH_LIBMESH
#include "MainOld.hpp"
#include "utopia_libmesh_Library.hpp"
#endif

#ifdef UTOPIA_WITH_MARS
#include "utopia_mars_Library.hpp"
#endif

int main(const int argc, char *argv[]) {
#ifdef UTOPIA_WITH_LIBMESH
    if (argc > 1 && argv[1] == std::string("--old")) {
        return MainOld(argc, argv);
    } else {
        utopia::Utopia::instance().add_library(utopia::make_unique<utopia::LibMeshLibrary>());
        return UTOPIA_MAIN(argc, argv);
    }
#endif  // UTOPIA_WITH_LIBMESH

    // #ifdef UTOPIA_WITH_MARS
    //     utopia::Utopia::instance().add_library(utopia::make_unique<utopia::MarsLibrary>());
    // #endif  // UTOPIA_WITH_MARS
}
