#include "utopia_stk_Mesh.hpp"

#include "utopia_make_unique.hpp"

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>

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

                // FIXME
                n_local_elements = n_elements;
                n_local_nodes = n_nodes;

                // const stk::topology theTopo = part->topology();
            }
        };

        Mesh::~Mesh() = default;
        Mesh::Mesh(const Comm &comm) : impl_(utopia::make_unique<Impl>()) { impl_->comm = comm; }

        bool Mesh::read(const Path &path) {
            auto comm = impl_->comm.raw_comm();

            auto meta_data = std::make_shared<Impl::MetaData>();
            auto bulk_data = std::make_shared<Impl::BulkData>(*meta_data, comm, Impl::BulkData::NO_AUTO_AURA);

            try {
                Impl::IOBroker io_broker(comm);
                io_broker.set_bulk_data(*bulk_data);
                // io_broker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "rcb"));
                io_broker.add_mesh_database(path.to_string(), ::stk::io::READ_MESH);
                io_broker.create_input_mesh();
                io_broker.populate_mesh();
                io_broker.populate_field_data();

                impl_->meta_data = meta_data;
                impl_->bulk_data = bulk_data;

                impl_->compute_mesh_stats();

                return true;
            } catch (const std::exception &ex) {
                utopia::err() << "Mesh::read(\"" << path.to_string() << "\") error: " << ex.what() << '\n';
                assert(false);
                return false;
            }
        }

        bool Mesh::write(const Path &path) const {
            auto comm = impl_->comm.raw_comm();
            auto &bulk_data = *impl_->bulk_data;

            try {
                Impl::IOBroker io_broker(comm);
                io_broker.set_bulk_data(bulk_data);
                auto out_id = io_broker.create_output_mesh(path.to_string(), ::stk::io::WRITE_RESTART);
                io_broker.write_output_mesh(out_id);
                return true;

            } catch (const std::exception &ex) {
                utopia::err() << "Mesh::write(\"" << path.to_string() << "\") error: " << ex.what() << '\n';
                return false;
            }
        }

        void Mesh::read(Input &in) {
            Path path;
            in.get("path", path);
            read(path);
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

    }  // namespace stk
}  // namespace utopia
