#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "utopia_libmesh_Library.hpp"

using namespace std;

int main(const int argc, char *argv[]) {
    // TODO
    // #ifdef UTOPIA_WITH_LIBMESH
    utopia::Utopia::instance().add_library(utopia::make_unique<utopia::LibMeshLibrary>());
    // #endif //UTOPIA_WITH_LIBMESH

    return UTOPIA_TEST(argc, argc);
}
