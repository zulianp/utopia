#include "utopia_Main.hpp"
#include "utopia_fe_config.hpp"

#ifdef UTOPIA_WITH_LIBMESH
#include "MainOld.hpp"
#include "utopia_libmesh_Library.hpp"
#endif

int main(const int argc, char *argv[]) {
#ifdef UTOPIA_WITH_LIBMESH
    if (argc > 1 && argv[1] == std::string("--old")) {
        return MainOld(argc, argv);
    } else {
        utopia::Utopia::instance().add_library(utopia::make_unique<utopia::LibMeshLibrary>());
    }
#endif  // UTOPIA_WITH_LIBMESH

    return UTOPIA_MAIN(argc, argv);
}
