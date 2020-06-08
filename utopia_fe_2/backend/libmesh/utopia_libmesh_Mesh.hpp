#ifndef UTOPIA_LIBMESH_MESH_HPP
#define UTOPIA_LIBMESH_MESH_HPP

// Utopia
#include "utopia_Communicator.hpp"
#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"

// FIXME
#include "utopia_petsc_Vector.hpp"

// Utopia FE
#include "utopia_Mesh.hpp"
#include "utopia_libmesh_ForwardDeclarations.hpp"

// Libmesh
#include "libmesh/libmesh_common.h"

namespace utopia {

    template <>
    class Mesh<libMesh::UnstructuredMesh> : public Configurable, public Describable {
    public:
        using Scalar = libMesh::Real;

        Mesh(const Communicator &comm);
        ~Mesh();

        void read(Input &is) override;
        void describe(std::ostream &os = std::cout) const override;

    private:
        std::shared_ptr<libMesh::Parallel::Communicator> comm_;
        std::unique_ptr<libMesh::UnstructuredMesh> impl_;
    };

    using LMMesh = utopia::Mesh<libMesh::UnstructuredMesh>;

    template <>
    class Traits<LMMesh> : public Traits<PetscVector> {};
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_MESH_HPP
