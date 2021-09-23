#ifndef UTOPIA_LIBMESH_MESH_HPP
#define UTOPIA_LIBMESH_MESH_HPP

#include "utopia_Communicator.hpp"
#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"

#include "utopia_fe_base.hpp"

#include "utopia_libmesh_ForwardDeclarations.hpp"

#include "utopia_libmesh_SideSet.hpp"

#include "utopia_AABB.hpp"

namespace utopia {

    template <>
    class Traits<utopia::libmesh::Mesh> : public Traits<UVector> {
    public:
        using Super = Traits<UVector>;
        using SideSet = utopia::libmesh::SideSet;
        using AABB = utopia::AABB<std::vector<Super::Scalar>>;
    };

    namespace libmesh {

        class Mesh final : public Configurable, public Describable {
        public:
            using Scalar = Traits<Mesh>::Scalar;
            using SizeType = Traits<Mesh>::SizeType;
            using Vector = Traits<Mesh>::Vector;
            using Comm = Traits<Mesh>::Communicator;
            using AABB = Traits<Mesh>::AABB;

            ~Mesh();
            Mesh(const Comm &comm = Comm::get_default());

            bool read(const Path &path);
            bool write(const Path &path);

            void read(Input &in) override;
            void describe(std::ostream &os) const override;

            const Comm &comm() const;
            Comm &comm();

            libMesh::MeshBase &raw_type();
            const libMesh::MeshBase &raw_type() const;
            void wrap(const std::shared_ptr<libMesh::MeshBase> &mesh);
            bool empty() const;

            void unit_cube(const SizeType &nx = 10, const SizeType &ny = 10, const SizeType &nz = 10);

            void box(const AABB &box,
                     const SizeType &nx,
                     const SizeType &ny,
                     const SizeType &nz,
                     const std::string &elem_type = "HEX8");

            int manifold_dimension() const;
            int spatial_dimension() const;

            SizeType n_nodes() const;
            SizeType n_local_nodes() const;

            void set_database(const Path &path);
            const Path &database() const;

            void displace(const Vector &displacement);

            void uniform_refine(const int n_refinements);

            void scale(const Scalar &scale_factor);

            SizeType n_elements() const;

            void bounding_box(AABB &output) const;

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;

            friend class MeshInitializer;
            friend class FunctionSpace;

            void init_distributed();
            void init_serial();
            void init_replicated();
        };

    }  // namespace libmesh

    using LMMesh = utopia::libmesh::Mesh;

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_MESH_HPP