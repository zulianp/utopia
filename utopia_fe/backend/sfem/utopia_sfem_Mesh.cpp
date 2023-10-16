#include "utopia_sfem_Mesh.hpp"

#include "utopia_IOStream.hpp"
#include "utopia_make_unique.hpp"

#include "read_mesh.h"
#include "sfem_base.h"

namespace utopia {

    namespace sfem {

        class Mesh::Impl {
        public:
            Communicator comm;
            mesh_t mesh;
            Impl() : comm(Communicator::get_default()) {}
        };

        Mesh::Mesh() : impl_(utopia::make_unique<Impl>()) {}
        Mesh::~Mesh() {}

        const Mesh::Communicator &Mesh::comm() const { return impl_->comm; }
        Mesh::Communicator &Mesh::comm() { return impl_->comm; }

        bool Mesh::read(const Path &path) {
            return SFEM_OK == mesh_read(impl_->comm.get(), path.c_str(), &impl_->mesh);
        }

        void Mesh::read(Input &in) {
            Path path;
            in.require("path", path);
            if (!this->read(path)) {
                utopia::err() << "Unable to read mesh at " << path << "\n";
            }
        }

        void Mesh::describe(std::ostream &os) const {
            // TODO
        }

    }  // namespace sfem
}  // namespace utopia