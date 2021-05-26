#include "utopia_mars_Mesh.hpp"

#include "utopia_make_unique.hpp"

#include "mars.hpp"
#include "mars_base.hpp"
#include "mars_context.hpp"
#include "mars_execution_context.hpp"

#include <functional>

namespace utopia {

    namespace mars {

        class Mesh::Impl {
        public:
            Comm comm;

            ::mars::proc_allocation resources;
            ::mars::context context;

            std::shared_ptr<::mars::DimensionlessIMesh> mesh;

            Impl(const Comm &comm) : comm(comm), context(::mars::make_context(resources, comm.raw_comm())) {}

            template <class MeshD>
            void wrap(const std::shared_ptr<MeshD> &mesh_d) {
                mesh = mesh_d;

                describe = [mesh_d]() { mesh_d->print_sfc(); };
            }

            std::function<int()> spatial_dimension;
            std::function<int()> manifold_dimension;
            std::function<SizeType()> n_elements;
            std::function<SizeType()> n_nodes;
            std::function<SizeType()> n_local_nodes;
            std::function<bool(const Path &path)> write;

            std::function<void()> describe;
        };

        Mesh::~Mesh() {}

        Mesh::Mesh(const Comm &comm) : impl_(utopia::make_unique<Impl>(comm)) {}

        void Mesh::read(Input &in) {
            std::string type = "cube";
            in.get("type", type);

            SizeType nx = 10;
            SizeType ny = 10;
            SizeType nz = 10;

            in.get("nx", nx);
            in.get("ny", ny);
            in.get("nz", nz);

            if (type == "square") {
                nz = 0;
            }

            unit_cube(nx, ny, nz);
        }

        void Mesh::describe(std::ostream &os) const {
            if (impl_->describe) {
                impl_->describe();
            } else {
                os << "utopia::mars::Mesh == null\n";
            }
        }

        const Mesh::Comm &Mesh::comm() const {}

        bool Mesh::empty() const { return !static_cast<bool>(impl_->mesh); }

        int Mesh::spatial_dimension() const {}

        Mesh::SizeType Mesh::n_elements() const {}

        Mesh::SizeType Mesh::n_nodes() const {}

        Mesh::SizeType Mesh::n_local_elements() const {}

        Mesh::SizeType Mesh::n_local_nodes() const {}

        void Mesh::unit_cube(const SizeType &nx, const SizeType &ny, const SizeType &nz) {
            // 2D
            if (nz == 0) {
                using DMesh = ::mars::DistributedMesh<::mars::ElementType::Quad4>;
                auto mesh = std::make_shared<DMesh>();
                ::mars::generate_distributed_cube(impl_->context, *mesh, nx, ny, nz);

                impl_->wrap(mesh);

            } else {
                // 3D
                assert(nx != 0);
                assert(ny != 0);
                assert(nz != 0);
                using DMesh = ::mars::DistributedMesh<::mars::ElementType::Hex8>;
                auto mesh = std::make_shared<DMesh>();
                ::mars::generate_distributed_cube(impl_->context, *mesh, nx, ny, nz);

                impl_->wrap(mesh);
            }
        }

        void Mesh::init() {}

    }  // namespace mars
}  // namespace utopia
