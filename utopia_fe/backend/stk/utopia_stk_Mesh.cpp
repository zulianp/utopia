#include "utopia_stk_Mesh.hpp"

#include "utopia_make_unique.hpp"

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

//  stk::ParallelMachine pm = FluyaEnv::self().parallel_comm();

//     metaData_ = std::make_unique<MetaData>();
//     bulkData_ = std::make_unique<BulkData>(*metaData_, pm, BulkData::NO_AUTO_AURA);
// //    bulkData_ = std::make_unique<BulkData>(*metaData_, pm, BulkData::AUTO_AURA);
//     ioBroker_ = std::make_unique<stk::io::StkMeshIoBroker>(pm);
//     ioBroker_->set_bulk_data(*bulkData_);

//     // untested, auto decomposition of the mesh
// //    ioBroker_->property_add(Ioss::Property("DECOMPOSITION_METHOD", "rcb"));

//     inputMeshIdx_ = ioBroker_->add_mesh_database(inputDBName_, stk::io::READ_MESH);
//     ioBroker_->create_input_mesh();

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
                io_broker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "rcb"));
                io_broker.add_mesh_database(path.to_string(), ::stk::io::READ_MESH);
                io_broker.create_input_mesh();
                io_broker.populate_mesh();
                io_broker.populate_field_data();

                impl_->meta_data = meta_data;
                impl_->bulk_data = bulk_data;

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
                // io_broker.set_active_mesh(out_id);
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

        void Mesh::describe(std::ostream &os) const {}

        const Mesh::Comm &Mesh::comm() const { return impl_->comm; }

        ::stk::mesh::BulkData &Mesh::raw_type() {
            assert(impl_->bulk_data);
            return *impl_->bulk_data;
        }

        const ::stk::mesh::BulkData &Mesh::raw_type() const {
            assert(impl_->bulk_data);
            return *impl_->bulk_data;
        }

        void Mesh::wrap(const std::shared_ptr<::stk::mesh::MetaData> &meta_data,
                        const std::shared_ptr<::stk::mesh::BulkData> &bulk_data) {
            impl_->meta_data = meta_data;
            impl_->bulk_data = bulk_data;
        }

        bool Mesh::empty() const { return !static_cast<bool>(impl_->bulk_data); }

    }  // namespace stk
}  // namespace utopia
