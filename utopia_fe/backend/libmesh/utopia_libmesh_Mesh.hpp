#ifndef UTOPIA_LIBMESH_MESH_HPP
#define UTOPIA_LIBMESH_MESH_HPP

#include "utopia_Communicator.hpp"
#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_Mesh.hpp"

namespace libMesh {
    class MeshBase;

    namespace Parallel {
        class Communicator;
    }

}  // namespace libMesh

namespace utopia {

    namespace libmesh {

        class MeshInitializer;

        class Mesh final : public Configurable, public Describable {
        public:
            Mesh();
            ~Mesh();
            Mesh(const Communicator &comm);

            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            libMesh::MeshBase &raw_type();
            const libMesh::MeshBase &raw_type() const;
            void wrap(const std::shared_ptr<libMesh::MeshBase> &mesh);
            bool empty() const;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;

            friend class MeshInitializer;
            void init_distributed();
            void init_serial();
            void init_replicated();
        };

    }  // namespace libmesh

    using LMMesh = utopia::libmesh::Mesh;

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_MESH_HPP