#ifndef UTOPIA_LIBMESH_MESH_HPP
#define UTOPIA_LIBMESH_MESH_HPP

#include "utopia_Communicator.hpp"
#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_Mesh.hpp"

#include "utopia_fe_base.hpp"

#include "utopia_libmesh_ForwardDeclarations.hpp"

namespace utopia {

    template <>
    class Traits<utopia::libmesh::Mesh> : public Traits<UVector> {};

    namespace libmesh {

        class Mesh final : public Configurable, public Describable {
        public:
            using SizeType = Traits<Mesh>::SizeType;
            using Comm = Traits<Mesh>::Communicator;

            ~Mesh();
            Mesh(const Comm &comm = Comm::get_default());

            bool read(const Path &path);
            bool write(const Path &path);

            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            const Comm &comm() const;

            libMesh::MeshBase &raw_type();
            const libMesh::MeshBase &raw_type() const;
            void wrap(const std::shared_ptr<libMesh::MeshBase> &mesh);
            bool empty() const;

            void unit_cube(const SizeType &nx = 10, const SizeType &ny = 10, const SizeType &nz = 10);

            int manifold_dimension() const;
            int spatial_dimension() const;

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