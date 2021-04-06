#include "utopia_stk_IO.hpp"

#include "utopia_stk_Mesh.hpp"

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace utopia {
    namespace stk {

        class IO::Impl : public Configurable {
        public:
            using MetaData = ::stk::mesh::MetaData;
            using BulkData = ::stk::mesh::BulkData;
            using IOBroker = ::stk::io::StkMeshIoBroker;

            void read(Input &in) override {
                /////////////////////////////////////////////
                bool decomposition_in_subdirs = false;
                in.get("decomposition_in_subdirs", decomposition_in_subdirs);
                in.get("decomposition_method", decomposition_method);
                in.get("verbose", verbose);
                in.get("import_all_field_data", import_all_field_data);

                // if (decomposition_in_subdirs) {
                //     decomposition_method = "";
                // }

                /////////////////////////////////////////////
                Path path;
                in.get("path", path);

                auto &comm = mesh.comm();

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

            bool load() {
                auto meta_data = std::make_shared<Impl::MetaData>();
                auto bulk_data = std::make_shared<Impl::BulkData>(*meta_data, mesh.comm().raw_comm(), auto_aura_option);

                try {
                    io_broker->set_bulk_data(*bulk_data);

                    if (!decomposition_method.empty()) {
                        io_broker->property_add(::Ioss::Property("DECOMPOSITION_METHOD", decomposition_method));
                    }

                    io_broker->add_mesh_database(read_specification, read_purpose);
                    io_broker->create_input_mesh();

                    if (import_all_field_data) {
                        io_broker->add_all_mesh_fields_as_input_fields();
                    }

                    io_broker->populate_mesh();

                    io_broker->populate_field_data();

                    mesh.wrap(meta_data, bulk_data);

                    mesh.init();

                    // if (verbose) {
                    //     std::vector<std::string> names;
                    //     io_broker->get_global_variable_names(names);

                    //     utopia::out() << "Variables:\n";
                    //     for (auto &n : names) {
                    //         utopia::out() << n << '\n';
                    //     }

                    //     utopia::out() << '\n';
                    // }

                    return true;
                } catch (const std::exception &ex) {
                    utopia::err() << "Mesh::read(\"" << read_specification << "\") error: " << ex.what() << '\n';
                    assert(false);
                    return false;
                }
            }

            bool write(const Path &write_path) {
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

            Impl(Mesh &mesh) : mesh(mesh), io_broker(utopia::make_unique<IOBroker>(mesh.comm().raw_comm())) {}

        public:
            Mesh &mesh;
            std::unique_ptr<IOBroker> io_broker;
            std::string read_specification;
            ::stk::io::DatabasePurpose read_purpose{::stk::io::READ_MESH};           //{READ_RESTART}
            BulkData::AutomaticAuraOption auto_aura_option{BulkData::NO_AUTO_AURA};  //{AUTO_AURA}
            std::string decomposition_method;
            bool create_edges{false};
            bool verbose{false};
            bool import_all_field_data{false};
        };

        void IO::read(Input &in) { impl_->read(in); }
        bool IO::load() { return impl_->load(); }
        bool IO::write(const Path &write_path) { return impl_->write(write_path); }

        void IO::set_read_path(const Path &path) { impl_->set_read_path(path); }
        void IO::set_read_specification(const std::string &format) { impl_->set_read_specification(format); }

        IO::IO(Mesh &mesh) : impl_(utopia::make_unique<Impl>(mesh)) {}

        IO::~IO() = default;

        void IO::import_all_field_data(const bool value) { impl_->import_all_field_data = value; }

    }  // namespace stk
}  // namespace utopia
