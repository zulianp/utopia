#include "utopia_stk_MeshIO.hpp"

#include "utopia_stk_Mesh.hpp"

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include "Trilinos_version.h"

#if TRILINOS_MAJOR_MINOR_VERSION >= 140100
#include <stk_mesh/base/MeshBuilder.hpp>
#endif

namespace utopia {
    namespace stk {

        class MeshIO::Impl : public Configurable {
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
                in.get("time_step_index", time_step_index);
                in.get("time", time_);
                in.get("output_path", output_path);

                Scalar scale = 1.;
                in.get("scale", scale);

                std::fill(scales.begin(), scales.end(), scale);
                in.get("scale_x", scales[0]);
                in.get("scale_y", scales[1]);
                in.get("scale_z", scales[2]);

                if (import_all_field_data) {
                    // read_purpose = ::stk::io::READ_RESTART;
                }

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

                    if (verbose) {
                        utopia::out() << "Reading DB at path: " << path << "\n";
                    }

                } else {
                    int nx = 10;
                    int ny = 10;
                    int nz = 10;

                    in.get("nx", nx);
                    in.get("ny", ny);
                    in.get("nz", nz);

                    Scalar box_min[3] = {0.0, 0.0, 0.0};
                    Scalar box_max[3] = {1.0, 1.0, 1.0};

                    in.get("min_x", box_min[0]);
                    in.get("min_y", box_min[1]);
                    in.get("min_z", box_min[2]);

                    in.get("max_x", box_max[0]);
                    in.get("max_y", box_max[1]);
                    in.get("max_z", box_max[2]);

                    std::string elem_type = "HEX8";
                    in.get("elem_type", elem_type);

                    if (verbose) {
                        utopia::out() << "Generating cube with " << nx << " x " << ny << " x " << nz << "\n";
                    }

                    // See packages/seacas/libraries/ioss/src/generated/Iogn_GeneratedMesh.C

                    std::string sidesets = "|sideset:xyzXYZ";
                    // std::string sidesets = "sideset:zyXYxZ";
                    std::string bbox = "bbox:";
                    bbox += std::to_string(box_min[0]) + ", ";
                    bbox += std::to_string(box_min[1]) + ", ";
                    bbox += std::to_string(box_min[2]) + ", ";
                    bbox += std::to_string(box_max[0]) + ", ";
                    bbox += std::to_string(box_max[1]) + ", ";
                    bbox += std::to_string(box_max[2]);

                    std::string specification = "generated:" + std::to_string(nx) + "x" + std::to_string(ny) + "x" +
                                                std::to_string(nz) + "|" + sidesets + "|" + bbox;

                    if (elem_type == "TET4") {
                        specification += "|tets";
                    }

                    in.get("specification", specification);
                    set_read_specification(specification);

                    mesh.set_is_generated_cube(true);
                }

                /////////////////////////////////////////////
                std::string read_purpose_str;
                in.get("read_purpose", read_purpose_str);
                if (read_purpose_str == "READ_RESTART") {
                    read_purpose = ::stk::io::READ_RESTART;
                }
                /////////////////////////////////////////////
                bool auto_aura = true;
                in.get("auto_aura", auto_aura);
                if (auto_aura) {
                    auto_aura_option = BulkData::AUTO_AURA;
                }
                /////////////////////////////////////////////
            }

            bool load_time_step(const Scalar t) {
                try {
                    io_broker->read_defined_input_fields(t);
                    return true;
                } catch (const std::exception &ex) {
                    utopia::err() << "Mesh::load_time_step(\"" << t << "\") error: " << ex.what() << '\n';
                    assert(false);
                    return false;
                }
            }

            bool load() {
                auto meta_data = std::make_shared<Impl::MetaData>();

#if TRILINOS_MAJOR_MINOR_VERSION >= 140100
                ::stk::mesh::MeshBuilder builder(mesh.comm().raw_comm());
                builder.set_aura_option(auto_aura_option);
                std::shared_ptr<Impl::BulkData> bulk_data = builder.create(meta_data);
#else
                auto bulk_data = std::make_shared<Impl::BulkData>(*meta_data, mesh.comm().raw_comm(), auto_aura_option);
#endif

                try {
                    io_broker->set_bulk_data(*bulk_data);

                    if (!decomposition_method.empty()) {
                        io_broker->property_add(::Ioss::Property("DECOMPOSITION_METHOD", decomposition_method));
                    }

                    io_broker->add_mesh_database(read_specification, read_purpose);
                    io_broker->create_input_mesh();

                    if (import_all_field_data) {
                        io_broker->add_all_mesh_fields_as_input_fields(time_match);
                    }

                    io_broker->populate_mesh();

                    io_broker->populate_field_data();

                    if (import_all_field_data) {
                        if (time_ < 0) {
                            io_broker->read_defined_input_fields(time_step_index);
                        } else {
                            io_broker->read_defined_input_fields(time_);
                        }
                    }

                    mesh.wrap(meta_data, bulk_data);

                    mesh.init();

                    input_id = io_broker->get_active_mesh();

                    if (scales[0] != 1. || scales[1] != 1. || scales[2] != 1.) {
                        mesh.scale(scales);
                    }

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
                    output_path = write_path;
                    io_broker->set_bulk_data(bulk_data);
                    output_id = io_broker->create_output_mesh(write_path.to_string(), write_purpose);
                    io_broker->write_output_mesh(output_id);
                    return true;

                } catch (const std::exception &ex) {
                    utopia::err() << "Mesh::write(\"" << write_path.to_string() << "\") error: " << ex.what() << '\n';
                    return false;
                }
            }

            void set_read_path(const Path &path) { this->read_specification = path.to_string(); }

            void set_read_specification(const std::string &format) { read_specification = format; }

            Impl(Mesh &mesh) : mesh(mesh), io_broker(utopia::make_unique<IOBroker>(mesh.comm().raw_comm())) {}

            inline bool empty_input() const { return input_id == -1; }

        public:
            Mesh &mesh;
            std::unique_ptr<IOBroker> io_broker;
            std::string read_specification;
            ::stk::io::DatabasePurpose read_purpose{::stk::io::READ_MESH};  //{READ_RESTART}
            ::stk::io::DatabasePurpose write_purpose{::stk::io::WRITE_RESTART};
            BulkData::AutomaticAuraOption auto_aura_option{BulkData::NO_AUTO_AURA};  //{AUTO_AURA}
            std::string decomposition_method;
            bool verbose{false};
            bool import_all_field_data{false};
            int time_step_index{1};
            double time_{-1};
            Path output_path{"./out.e"};
            int output_id{-1};
            int input_id{-1};
            ::stk::io::MeshField::TimeMatchOption time_match{::stk::io::MeshField::CLOSEST};
            std::vector<Scalar> scales{{1, 1, 1}};
        };

        void MeshIO::enable_interpolation_mode() { impl_->time_match = ::stk::io::MeshField::LINEAR_INTERPOLATION; }

        void MeshIO::read(Input &in) { impl_->read(in); }
        bool MeshIO::load() { return impl_->load(); }
        bool MeshIO::write(const Path &write_path) { return impl_->write(write_path); }

        void MeshIO::set_read_path(const Path &path) { impl_->set_read_path(path); }
        void MeshIO::set_read_specification(const std::string &format) { impl_->set_read_specification(format); }

        MeshIO::MeshIO(Mesh &mesh) : impl_(utopia::make_unique<Impl>(mesh)) {}

        MeshIO::~MeshIO() = default;

        void MeshIO::import_all_field_data(const bool value) { impl_->import_all_field_data = value; }

        void MeshIO::set_output_path(const Path &path) { impl_->output_path = path; }

        void MeshIO::set_output_mode(enum OutputMode output_mode) {
            switch (output_mode) {
                case OUTPUT_MODE_OVERWRITE: {
                    impl_->write_purpose = ::stk::io::WRITE_RESULTS;
                    break;
                }
                case OUTPUT_MODE_APPEND: {
                    impl_->write_purpose = ::stk::io::APPEND_RESULTS;
                    break;
                }
                case OUTPUT_MODE_RESTART: {
                    impl_->write_purpose = ::stk::io::WRITE_RESTART;
                    break;
                }

                default:
                    Utopia::Abort("Unsupported output mode");
            }
        }

        void MeshIO::set_output_mode(const std::string &output_mode) {
            if (output_mode.empty()) return;

            if (output_mode == "APPEND") {
                set_output_mode(OUTPUT_MODE_APPEND);
            } else if (output_mode == "OVERWRITE") {
                set_output_mode(OUTPUT_MODE_OVERWRITE);
            } else if (output_mode == "RESTART") {
                set_output_mode(OUTPUT_MODE_RESTART);
            } else {
                Utopia::Abort("Unsupported output mode: " + output_mode + "! Use APPEND or OVERWRITE");
            }
        }

        const Path &MeshIO::output_path() const { return impl_->output_path; }

        ::stk::io::StkMeshIoBroker &MeshIO::raw_type() {
            assert(impl_->io_broker);
            return *impl_->io_broker;
        }

        void MeshIO::create_output_mesh() {
            assert(impl_->output_id == -1);
            impl_->io_broker->set_bulk_data(impl_->mesh.bulk_data());

            impl_->output_id =
                impl_->io_broker->create_output_mesh(impl_->output_path.to_string(), impl_->write_purpose);
        }

        bool MeshIO::ensure_output() {
            if (impl_->output_id != -1) {
                return true;
            }

            create_output_mesh();
            return true;
        }

        void MeshIO::register_output_field(const std::string &var_name) {
            if (ensure_output()) {
                assert(impl_->output_id != -1);
                ::stk::mesh::FieldBase *field = ::stk::mesh::get_field_by_name(var_name, impl_->mesh.meta_data());
                if (field) {
                    impl_->io_broker->add_field(impl_->output_id, *field, var_name);
                } else {
                    utopia::err() << "[Error] mesh does not have a field with name " << var_name << '\n';
                    assert(false);
                }
            }
        }

        bool MeshIO::write(const int step, const Scalar t) {
            UTOPIA_UNUSED(step);

            try {
                assert(impl_->output_id != -1);
                impl_->io_broker->process_output_request(impl_->output_id, t);
                return true;
            } catch (const std::exception &ex) {
                utopia::err() << "MeshIO::write: " << impl_->output_path.to_string() << " error: " << ex.what() << '\n';
                assert(false);
                return false;
            }
        }

        std::vector<double> MeshIO::get_time_steps() const { return impl_->io_broker->get_time_steps(); }

        bool MeshIO::load_time_step(const Scalar t) { return impl_->load_time_step(t); }

        int MeshIO::num_time_steps() const { return impl_->io_broker->get_num_time_steps(); }
        MeshIO::Scalar MeshIO::max_time() const { return impl_->io_broker->get_max_time(); }

        bool MeshIO::load_last_time_step() { return load_time_step(max_time()); }

        void MeshIO::set_import_all_data(const bool val) { impl_->import_all_field_data = val; }

        bool MeshIO::ensure_input() {
            if (impl_->empty_input()) {
                return impl_->load();
            }

            return true;
        }

    }  // namespace stk
}  // namespace utopia
