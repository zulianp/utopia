#include "utopia_stk_Mesh.hpp"

#include "utopia_make_unique.hpp"

#include "utopia_stk_Commons.hpp"

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

            friend class Mesh::IO;
        };

        class Mesh::IO : public Configurable {
        public:
            using MetaData = ::stk::mesh::MetaData;
            using BulkData = ::stk::mesh::BulkData;
            using IOBroker = ::stk::io::StkMeshIoBroker;

            void read(Input &in) override {
                /////////////////////////////////////////////
                bool decomposition_in_subdirs = false;
                in.get("decomposition_in_subdirs", decomposition_in_subdirs);
                in.get("decomposition_method", decomposition_method);

                if (!decomposition_method.empty()) {
                    assert(!decomposition_in_subdirs && "Either one or the other!");
                }

                /////////////////////////////////////////////
                Path path;
                in.get("path", path);

                if (!path.empty()) {
                    if (decomposition_in_subdirs && comm.size() > 1) {
                        read_specification = path.parent() / Path(std::to_string(comm.size())) /
                                             (path.file_name() + "." + path.extension());
                    } else {
                        read_specification = path;
                    }
                } else {
                    std::string specification = "generated:10x10x10|sidesetxX";
                    in.get("specification", specification);
                    set_read_specification(specification);
                }

                /////////////////////////////////////////////
                std::string read_purpose_str;
                in.get("read_purpose", read_purpose_str);
                if (read_purpose_str == "READ_RESTART") {
                    read_purpose = ::stk::io::READ_RESTART;
                }
                /////////////////////////////////////////////
                bool auto_aura = false;
                in.get("auto_aura", auto_aura);
                if (auto_aura) {
                    auto_aura_option = BulkData::AUTO_AURA;
                }
                /////////////////////////////////////////////
            }

            bool load(Mesh &mesh) {
                auto meta_data = std::make_shared<Impl::MetaData>();
                auto bulk_data = std::make_shared<Impl::BulkData>(*meta_data, mesh.comm().raw_comm(), auto_aura_option);

                try {
                    io_broker->set_bulk_data(*bulk_data);
                    if (!decomposition_method.empty()) {
                        io_broker->property_add(::Ioss::Property("DECOMPOSITION_METHOD", decomposition_method));
                    }

                    io_broker->add_mesh_database(read_specification, read_purpose);
                    io_broker->create_input_mesh();
                    io_broker->populate_mesh();
                    io_broker->populate_field_data();

                    mesh.wrap(meta_data, bulk_data);

                    mesh.impl_->compute_mesh_stats();

                    return true;
                } catch (const std::exception &ex) {
                    utopia::err() << "Mesh::read(\"" << read_specification << "\") error: " << ex.what() << '\n';
                    assert(false);
                    return false;
                }
            }

            bool write(const Path &write_path, const Mesh &mesh) {
                auto &bulk_data = mesh.bulk_data();

                try {
                    io_broker->set_bulk_data(bulk_data);
                    auto out_id = io_broker->create_output_mesh(write_path.to_string(), ::stk::io::WRITE_RESTART);
                    io_broker->write_output_mesh(out_id);
                    return true;

                } catch (const std::exception &ex) {
                    utopia::err() << "Mesh::write(\"" << write_path.to_string() << "\") error: " << ex.what() << '\n';
                    return false;
                }
            }

            void set_read_path(const Path &path) { this->read_specification = path.to_string(); }

            void set_read_specification(const std::string &format) { read_specification = format; }

            IO(const Communicator &comm) : comm(comm), io_broker(utopia::make_unique<IOBroker>(comm.raw_comm())) {}

        public:
            const Communicator &comm;
            std::unique_ptr<IOBroker> io_broker;
            std::string read_specification;
            ::stk::io::DatabasePurpose read_purpose{::stk::io::READ_MESH};           //{READ_RESTART}
            BulkData::AutomaticAuraOption auto_aura_option{BulkData::NO_AUTO_AURA};  //{AUTO_AURA}
            std::string decomposition_method;
            bool create_edges{false};
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
            IO io(comm());
            io.set_read_path(path);
            return io.load(*this);
        }

        bool Mesh::write(const Path &path) const {
            IO io(comm());
            return io.write(path, *this);
        }

        void Mesh::read(Input &in) {
            IO io(comm());
            io.read(in);
            if (!io.load(*this)) {
                assert(false);
            }
        }

        void Mesh::unit_cube(const SizeType &nx, const SizeType &ny, const SizeType &nz) {
            IO io(comm());

            char format[100];
            std::sprintf(format, "generated:%dx%dx%d|sideset:xX", int(nx), int(ny), int(nz));
            io.set_read_specification(format);
            if (!io.load(*this)) {
                assert(false);
            }
        }

        void Mesh::describe(std::ostream &os) const { impl_->bulk_data->dump_all_mesh_info(os); }

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

    }  // namespace stk
}  // namespace utopia
