#include "utopia_stk_Mesh.hpp"

#include "utopia_make_unique.hpp"

#include "utopia_stk_Commons.hpp"

#include "utopia_stk_MeshIO.hpp"

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <cstdio>

namespace utopia {
    namespace stk {

        class Mesh::Impl {
        public:
            using MetaData = ::stk::mesh::MetaData;
            using BulkData = ::stk::mesh::BulkData;
            using IOBroker = ::stk::io::StkMeshIoBroker;

            Comm comm;
            std::shared_ptr<MetaData> meta_data;
            std::shared_ptr<BulkData> bulk_data;

            SizeType n_elements{-1};
            SizeType n_nodes{-1};
            SizeType n_local_elements{-1};
            SizeType n_local_nodes{-1};

            void compute_mesh_stats() {
                std::vector<size_t> entity_counts;
                ::stk::mesh::comm_mesh_counts(*bulk_data, entity_counts);

                n_elements = entity_counts[::stk::topology::ELEMENT_RANK];
                n_nodes = entity_counts[::stk::topology::NODE_RANK];

                n_local_elements = utopia::stk::count_local_elements(*bulk_data);
                n_local_nodes = utopia::stk::count_local_nodes(*bulk_data);
            }
        };

        Mesh::~Mesh() = default;
        Mesh::Mesh(const Comm &comm) : impl_(utopia::make_unique<Impl>()) { impl_->comm = comm; }

        // 1736
        // https://github.com/NaluCFD/Nalu/blob/master/src/Realm.C

        // automatic_decomposition_type (DECOMPOSITION_METHOD)
        // Used only for parallel runs, this indicates how the a single mesh database must be decomposed amongst the MPI
        // processes during initialization. This option should not be used if the mesh has already been decomposed by an
        // external utility. Possible values are: Purpose Generate statistics for the flow field Extract integrated data
        // from the simulation Compare the solution error to a reference solution Extract data using probes Model
        // turbine blades/tower using actuator lines Momentum source term to drive ABL flows to a desired velocity
        // profile
        //         The name of the realm. The name provided here is used in the Time_Integrators section to determine
        //         the time-integration scheme used for this computational domain.
        // The name of the Exodus-II mesh file that defines the computational domain for this realm. Note that only the
        // base name (i.e., without the .NPROCS.IPROC suffix) is provided even for simulations using previously
        // decomposed mesh/restart files. Value Description rcb recursive coordinate bisection rib recursive inertial
        // bisection linear elements in order first n/p to proc 0, next to proc 1. cyclic elements handed out to id %
        // proc_count

        bool Mesh::read(const Path &path) {
            MeshIO io(*this);
            io.set_read_path(path);
            return io.load();
        }

        bool Mesh::write(const Path &path) {
            MeshIO io(*this);
            return io.write(path);
        }

        void Mesh::read(Input &in) {
            MeshIO io(*this);
            io.read(in);
            if (!io.load()) {
                assert(false);
            }
        }

        void Mesh::unit_cube(const SizeType &nx, const SizeType &ny, const SizeType &nz) {
            MeshIO io(*this);

            char format[100];
            std::sprintf(format, "generated:%dx%dx%d|sideset:xX", int(nx), int(ny), int(nz));
            io.set_read_specification(format);
            if (!io.load()) {
                assert(false);
            }
        }

        void Mesh::describe(std::ostream &os) const {
            if (comm().rank() == 0) {
                os << "Parts:\n";
                for (auto ptr : meta_data().get_parts()) {
                    auto &p = *ptr;
                    if (p.id() != -1) {
                        os << p.name() << ' ' << p.id() << '\n';
                    }
                }

                for (auto &field : impl_->meta_data->get_fields()) {
                    os << field->name() << ", " << field->entity_rank() << ", num states: " << field->number_of_states()
                       << '\n';
                }
            }

            // impl_->bulk_data->dump_all_mesh_info(os);
        }

        const Mesh::Comm &Mesh::comm() const { return impl_->comm; }

        ::stk::mesh::BulkData &Mesh::bulk_data() const {
            assert(impl_->bulk_data);
            return *impl_->bulk_data;
        }

        ::stk::mesh::MetaData &Mesh::meta_data() const {
            assert(impl_->meta_data);
            return *impl_->meta_data;
        }

        void Mesh::wrap(const std::shared_ptr<::stk::mesh::MetaData> &meta_data,
                        const std::shared_ptr<::stk::mesh::BulkData> &bulk_data) {
            impl_->meta_data = meta_data;
            impl_->bulk_data = bulk_data;

            if (!empty()) {
                impl_->compute_mesh_stats();
            }
        }

        bool Mesh::empty() const { return !static_cast<bool>(impl_->bulk_data); }

        int Mesh::spatial_dimension() const {
            assert(impl_->meta_data);
            return impl_->meta_data->spatial_dimension();
        }

        Mesh::SizeType Mesh::n_elements() const { return impl_->n_elements; }

        Mesh::SizeType Mesh::n_nodes() const { return impl_->n_nodes; }

        Mesh::SizeType Mesh::n_local_elements() const { return impl_->n_local_elements; }

        Mesh::SizeType Mesh::n_local_nodes() const { return impl_->n_local_nodes; }

        void Mesh::displace(const Vector &displacement) { assert(false); }

        void Mesh::init() { impl_->compute_mesh_stats(); }

    }  // namespace stk
}  // namespace utopia
