#include "utopia_moonolith_Mesh.hpp"

#include "utopia_Options.hpp"

// Moonolith includes
#include "moonolith_Mesh.hpp"
#include "moonolith_tri_mesh_reader.hpp"
#include "moonolith_vtu_mesh_writer.hpp"
#include "par_moonolith.hpp"

#include <functional>

namespace utopia {

    namespace moonolith {

        class Mesh::Impl {
        public:
            using IMesh = ::moonolith::IMesh<Scalar>;
            using Mesh1 = ::moonolith::Mesh<Mesh::Scalar, 1>;
            using Mesh2 = ::moonolith::Mesh<Mesh::Scalar, 2>;
            using Mesh3 = ::moonolith::Mesh<Mesh::Scalar, 3>;

            Comm comm;
            std::shared_ptr<IMesh> mesh;

            template <class MeshD>
            void wrap(const std::shared_ptr<MeshD> &mesh_d) {
                this->mesh = mesh_d;

                spatial_dimension = []() -> int { return MeshD::Dim; };
                manifold_dimension = [mesh_d]() -> int { return mesh_d->manifold_dim(); };
                n_elements = [mesh_d]() -> SizeType { return mesh_d->n_elements(); };
                n_nodes = [mesh_d, this]() -> SizeType {
                    // FIXME
                    return comm.sum(mesh_d->n_nodes());
                };

                n_local_nodes = [mesh_d]() -> SizeType {
                    // FIXME
                    return mesh_d->n_nodes();
                };

                write = [mesh_d](const Path &path) -> bool {
                    ::moonolith::VTUMeshWriter<MeshD> writer;
                    return writer.write(path.to_string(), *mesh_d);
                };
            }

            std::function<int()> spatial_dimension;
            std::function<int()> manifold_dimension;
            std::function<SizeType()> n_elements;
            std::function<SizeType()> n_nodes;
            std::function<SizeType()> n_local_nodes;
            std::function<bool(const Path &path)> write;
        };

        Mesh::~Mesh() = default;

        Mesh::Mesh(const Comm &comm) : impl_(utopia::make_unique<Impl>()) { impl_->comm = comm; }

        bool Mesh::read(const Path &path) {
            ::moonolith::TriMeshReader<Impl::Mesh3> reader;
            auto mesh = std::make_shared<Impl::Mesh3>(impl_->comm.raw_comm());

            if (reader.read(path.to_string(), *mesh)) {
                mesh->set_manifold_dim(2);
                impl_->wrap(mesh);
                return true;
            } else {
                return false;
            }
        }

        bool Mesh::write(const Path &path) const { return impl_->write(path); }

        void Mesh::read(Input &in) {
            Path path;
            if (!Options().add_option("path", path, "Path to the mesh. Format supported is .tri.").parse(in)) {
                return;
            }

            if (!read(path)) {
                return;
            }
        }

        void Mesh::describe(std::ostream &os) const {}

        const Mesh::Comm &Mesh::comm() const { return impl_->comm; }

        template <int Dim>
        const std::shared_ptr<::moonolith::Mesh<Mesh::Scalar, Dim>> Mesh::raw_type() const {
            if (Dim != impl_->spatial_dimension()) {
                assert(false && "Trying to read wrong dimension");
            }

            return std::dynamic_pointer_cast<::moonolith::Mesh<Mesh::Scalar, Dim>>(impl_->mesh);
        }

        template <int Dim>
        void Mesh::wrap(const std::shared_ptr<::moonolith::Mesh<Scalar, Dim>> &mesh) {
            impl_->wrap(mesh);
        }

        bool Mesh::empty() const { return !static_cast<bool>(impl_->mesh); }

        template void Mesh::wrap(const std::shared_ptr<::moonolith::Mesh<Mesh::Scalar, 1>> &);
        template void Mesh::wrap(const std::shared_ptr<::moonolith::Mesh<Mesh::Scalar, 2>> &);
        template void Mesh::wrap(const std::shared_ptr<::moonolith::Mesh<Mesh::Scalar, 3>> &);

        template const std::shared_ptr<::moonolith::Mesh<Mesh::Scalar, 1>> Mesh::raw_type() const;
        template const std::shared_ptr<::moonolith::Mesh<Mesh::Scalar, 2>> Mesh::raw_type() const;
        template const std::shared_ptr<::moonolith::Mesh<Mesh::Scalar, 3>> Mesh::raw_type() const;

        int Mesh::spatial_dimension() const {
            assert(impl_->mesh);
            return impl_->spatial_dimension();
        }

        int Mesh::manifold_dimension() const {
            assert(impl_->mesh);
            return impl_->manifold_dimension();
        }

        Mesh::SizeType Mesh::n_nodes() const { impl_->n_nodes(); }

        Mesh::SizeType Mesh::n_local_nodes() const { impl_->n_local_nodes(); }

    }  // namespace moonolith
}  // namespace utopia
