#ifndef UTOPIA_LIBMESH_LIBRARY_HPP
#define UTOPIA_LIBMESH_LIBRARY_HPP

#include "utopia_Library.hpp"

#include "utopia_libmesh_ForwardDeclarations.hpp"

namespace utopia {
    class LibMeshLibrary : public Library {
    public:
        void init(int argc, char *argv[]) override;
        int finalize() override;

        std::string name() const override;
        std::string version_string() const override;

        LibMeshLibrary();
        ~LibMeshLibrary() override;

    private:
        std::unique_ptr<libMesh::LibMeshInit> init_;
    };
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_LIBRARY_HPP