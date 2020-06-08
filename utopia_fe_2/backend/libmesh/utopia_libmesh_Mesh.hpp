#ifndef UTOPIA_LIBMESH_MESH_HPP
#define UTOPIA_LIBMESH_MESH_HPP

#include "utopia_Communicator.hpp"
#include "utopia_Input.hpp"

// FE stuff
#include "utopia_Mesh.hpp"
#include "utopia_libmesh_ForwardDeclarations.hpp"

namespace utopia {

    template <>
    class Mesh<libMesh::MeshBase> : public Configurable {
    public:
        Mesh(const Communicator &comm);
        ~Mesh();

        void read(Input &is) override;

    private:
        std::shared_ptr<libMesh::Parallel::Communicator> comm_;
        std::unique_ptr<libMesh::MeshBase> impl_;
    };

    using LMMesh = utopia::Mesh<libMesh::MeshBase>;
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_MESH_HPP
