#include "utopia_libmesh_Library.hpp"

#include "utopia_make_unique.hpp"

#include "utopia_fe_base.hpp"

#include "libmesh/libmesh.h"

namespace utopia {

    class LibMeshLibrary::Impl {
    public:
        Impl(int argc, char *argv[], MPI_Comm comm) : init(argc, argv, comm) {}

        libMesh::LibMeshInit init;
    };

    LibMeshLibrary::LibMeshLibrary() {}
    LibMeshLibrary::~LibMeshLibrary() {}
    void LibMeshLibrary::init(int argc, char *argv[]) {
        impl_ = utopia::make_unique<Impl>(argc, argv, Traits<UVector>::Communicator::get_default().raw_comm());
    }

    int LibMeshLibrary::finalize() {
        impl_.release();
        return 0;
    }
    std::string LibMeshLibrary::name() const { return "libMesh"; }
    std::string LibMeshLibrary::version_string() const { return LIBMESH_LIB_VERSION; }

}  // namespace utopia
