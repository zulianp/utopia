#include "utopia.hpp"
#include "utopia_Testing.hpp"

#include "utopia_fe_config.hpp"

#ifdef UTOPIA_WITH_LIBMESH
#include "utopia_libmesh_Library.hpp"
#endif  // UTOPIA_WITH_LIBMESH

int main(const int argc, char *argv[]) {
#ifdef UTOPIA_WITH_LIBMESH
    utopia::Utopia::instance().add_library(utopia::make_unique<utopia::LibMeshLibrary>());
#endif  // UTOPIA_WITH_LIBMESH

    return UTOPIA_TEST(argc, argc);
}
