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

            void init(const Communicator &in_comm) {
                comm_do_not_use_ = std::make_shared<libMesh::Parallel::Communicator>(in_comm.raw_comm());
            }

            inline bool empty() const { return !mesh; }

            std::shared_ptr<libMesh::Parallel::Communicator> comm_do_not_use_;
        };

        Mesh::~Mesh() {}

        Mesh::Mesh(const Communicator &comm) : impl_(utopia::make_unique<Impl>()) { impl_->init(comm); }

        void Mesh::read(Input &in) {
            UTOPIA_UNUSED(in);

            MeshInitializer init(*this);
            init.read(in);

            // If it is still empty initialize an empty mesh
            if (empty()) {
                init_distributed();
            }
        }

        bool Mesh::empty() const { return impl_->empty(); }

        void Mesh::describe(std::ostream &os) const { UTOPIA_UNUSED(os); }

        libMesh::MeshBase &Mesh::raw_type() { return *impl_->mesh; }

        const libMesh::MeshBase &Mesh::raw_type() const { return *impl_->mesh; }

        void Mesh::wrap(const std::shared_ptr<libMesh::MeshBase> &mesh) { impl_->mesh = mesh; }

        void Mesh::init_distributed() {
            impl_->mesh = std::make_shared<libMesh::DistributedMesh>(*impl_->comm_do_not_use_);
        }

        void Mesh::init_serial() { impl_->mesh = std::make_shared<libMesh::SerialMesh>(*impl_->comm_do_not_use_); }

        void Mesh::init_replicated() {
            impl_->mesh = std::make_shared<libMesh::ReplicatedMesh>(*impl_->comm_do_not_use_);
        }

        void Mesh::read(const Path &path) {
            if (empty()) {
                init_distributed();
            }

            impl_->mesh->read(path.to_string());
        }

        void Mesh::write(const Path &path) { libMesh::NameBasedIO(*impl_->mesh).write(path.to_string()); }

    }  // namespace libmesh

}  // namespace utopia
