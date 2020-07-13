#include "utopia_libmesh_Library.hpp"
#include "utopia_Base.hpp"

#include "libmesh/libmesh.h"
#include "petscsys.h"
#include "utopia_make_unique.hpp"

namespace utopia {

    LibMeshLibrary::LibMeshLibrary() {}

    LibMeshLibrary::~LibMeshLibrary() {}

    void LibMeshLibrary::init(int argc, char *argv[]) {
        init_ = utopia::make_unique<libMesh::LibMeshInit>(argc, argv, PETSC_COMM_WORLD);
    }

    int LibMeshLibrary::finalize() {
        init_ = nullptr;
        return 0;
    }

    std::string LibMeshLibrary::name() const { return "LibMesh"; }

    std::string LibMeshLibrary::version_string() const { return "TODO"; }
}  // namespace utopia
