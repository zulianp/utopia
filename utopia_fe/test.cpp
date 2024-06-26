#include "utopia.hpp"
#include "utopia_Testing.hpp"

#include "utopia_fe_config.hpp"

#ifdef UTOPIA_ENABLE_LIBMESH
#include "utopia_libmesh_Library.hpp"
#endif  // UTOPIA_ENABLE_LIBMESH

int main(const int argc, char *argv[]) {
#ifdef UTOPIA_ENABLE_LIBMESH
    utopia::Utopia::instance().add_library(utopia::make_unique<utopia::LibMeshLibrary>());
#endif  // UTOPIA_ENABLE_LIBMESH

    // #ifdef UTOPIA_ENABLE_MARS
    //     utopia::Utopia::instance().add_library(utopia::make_unique<utopia::MarsLibrary>());
    // #endif  // UTOPIA_ENABLE_MARS

    return UTOPIA_TEST(argc, argc);
}
