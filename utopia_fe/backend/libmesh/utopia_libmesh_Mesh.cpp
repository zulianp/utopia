#include "utopia_libmesh_Mesh.hpp"

#include "utopia_make_unique.hpp"

#include "utopia_fe_base.hpp"

#include "utopia_libmesh_MeshInitializer.hpp"

// All libmesh includes
#include "libmesh/mesh_base.h"
#include "libmesh/namebased_io.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/serial_mesh.h"

namespace utopia {

    namespace libmesh {
        class Mesh::Impl {
        public:
            std::shared_ptr<libMesh::MeshBase> mesh;

            void init(const Comm &in_comm) {
                comm = in_comm;
                comm_do_not_use_ = std::make_shared<libMesh::Parallel::Communicator>(in_comm.raw_comm());
            }

            inline bool empty() const { return !mesh; }

            Comm comm;
            std::shared_ptr<libMesh::Parallel::Communicator> comm_do_not_use_;
        };

        Mesh::~Mesh() {}

        Mesh::Mesh(const Comm &comm) : impl_(utopia::make_unique<Impl>()) { impl_->init(comm); }

        void Mesh::read(Input &in) {
            UTOPIA_UNUSED(in);

            MeshInitializer init(*this);
            init.read(in);

            // If it is still empty initialize an empty mesh
            if (empty()) {
                init_distributed();
            }
        }

        const Mesh::Comm &Mesh::comm() const { return impl_->comm; }

        bool Mesh::empty() const { return impl_->empty(); }

        void Mesh::describe(std::ostream &os) const { UTOPIA_UNUSED(os); }

        libMesh::MeshBase &Mesh::raw_type() {
            assert(impl_->mesh);
            return *impl_->mesh;
        }

        const libMesh::MeshBase &Mesh::raw_type() const {
            assert(impl_->mesh);
            return *impl_->mesh;
        }

        void Mesh::wrap(const std::shared_ptr<libMesh::MeshBase> &mesh) {
            impl_->comm.set(mesh->comm().get());
            impl_->mesh = mesh;
        }

        void Mesh::init_distributed() {
            assert(impl_->comm_do_not_use_);
            impl_->mesh = std::make_shared<libMesh::DistributedMesh>(*impl_->comm_do_not_use_);
        }

        void Mesh::init_serial() {
            assert(impl_->comm_do_not_use_);
            impl_->mesh = std::make_shared<libMesh::SerialMesh>(*impl_->comm_do_not_use_);
        }

        void Mesh::init_replicated() {
            assert(impl_->comm_do_not_use_);
            impl_->mesh = std::make_shared<libMesh::ReplicatedMesh>(*impl_->comm_do_not_use_);
        }

        void Mesh::read(const Path &path) {
            if (empty()) {
                init_distributed();
            }

            impl_->mesh->read(path.to_string());
        }

        void Mesh::write(const Path &path) { libMesh::NameBasedIO(*impl_->mesh).write(path.to_string()); }

        void Mesh::unit_cube(const SizeType &nx, const SizeType &ny, const SizeType &nz) {
            InputParameters params;
            params.set("type", "cube");
            params.set("nx", nx);
            params.set("ny", ny);
            params.set("nz", nz);

            MeshInitializer init(*this);
            init.read(params);
        }
    }  // namespace libmesh

}  // namespace utopia
