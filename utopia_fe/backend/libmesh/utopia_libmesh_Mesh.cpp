#include "utopia_libmesh_Mesh.hpp"

#include "utopia_make_unique.hpp"

#include "utopia_fe_base.hpp"

#include "utopia_libmesh_MeshInitializer.hpp"

// All libmesh includes
#include "libmesh/boundary_info.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/namebased_io.h"
#include "libmesh/node.h"
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
            Path database;
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

        void Mesh::describe(std::ostream &os) const {
            if (empty()) {
                os << "empty\n";
                return;
            }

            os << "Database: " << database() << '\n';

            auto &mesh = this->raw_type();
            auto &bi = mesh.get_boundary_info();

            os << "------------------\n";
            os << "Boundary info\n";
            os << "- Sidesets:\n";
            for (auto id : bi.get_side_boundary_ids()) {
                auto name = bi.get_sideset_name(id);
                os << '\t' << name << ": " << id << '\n';
            }

            if (!bi.get_edge_boundary_ids().empty()) {
                os << "- Edgesets:\n";

                for (auto id : bi.get_edge_boundary_ids()) {
                    auto name = bi.get_edgeset_name(id);
                    os << '\t' << name << ": " << id << '\n';
                }
            }

            if (!bi.get_node_boundary_ids().empty()) {
                os << "- Nodesets:\n";

                for (auto id : bi.get_node_boundary_ids()) {
                    auto name = bi.get_nodeset_name(id);
                    os << '\t' << name << ": " << id << '\n';
                }
            }

            os << "------------------\n";
        }

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

        bool Mesh::read(const Path &path) {
            if (empty()) {
                init_distributed();
            }

            impl_->mesh->read(path.to_string());
            set_database(path);
            return true;
        }

        bool Mesh::write(const Path &path) {
            // if (path.extension() == "e") {
            //     libMesh::ExodusII_IO(*impl_->mesh).write(path.to_string());
            // } else {
            libMesh::NameBasedIO(*impl_->mesh).write(path.to_string());
            // }
            return true;
        }

        void Mesh::unit_cube(const SizeType &nx, const SizeType &ny, const SizeType &nz) {
            InputParameters params;
            params.set("type", "cube");
            params.set("nx", nx);
            params.set("ny", ny);
            params.set("nz", nz);

            MeshInitializer init(*this);
            init.read(params);
        }

        int Mesh::manifold_dimension() const { return impl_->mesh->mesh_dimension(); }

        int Mesh::spatial_dimension() const { return impl_->mesh->spatial_dimension(); }

        Mesh::SizeType Mesh::n_nodes() const { return impl_->mesh->n_nodes(); }

        Mesh::SizeType Mesh::n_local_nodes() const { return impl_->mesh->n_local_nodes(); }

        Mesh::SizeType Mesh::n_elements() const { return impl_->mesh->n_active_elem(); }

        void Mesh::set_database(const Path &path) { impl_->database = path; }
        const Path &Mesh::database() const { return impl_->database; }

        void Mesh::displace(const Vector &displacement) {
            assert(displacement.comm().size() == 1);

            // FIXME
            std::vector<int> vars = {0, 1, 2};
            int sys_num = 0;

            Read<Vector> r_d(displacement);

            auto &mesh = raw_type();
            auto m_it = mesh.local_nodes_begin();
            auto m_end = mesh.local_nodes_end();

            const static int dim = spatial_dimension();

            for (; m_it != m_end; ++m_it) {
                for (int c = 0; c < dim; ++c) {
                    const auto dof_id = (*m_it)->dof_number(sys_num, vars[c], 0);
                    (**m_it)(c) += displacement.get(dof_id);
                }
            }
        }

        void Mesh::uniform_refine(const int n_refinements) {
            if (n_refinements <= 0) return;

            libMesh::MeshRefinement mesh_refinement(raw_type());
            mesh_refinement.make_flags_parallel_consistent();
            mesh_refinement.uniformly_refine(n_refinements);
        }

        void Mesh::scale(const Scalar &) { assert(false && "IMPLEMENT ME"); }

    }  // namespace libmesh

}  // namespace utopia
