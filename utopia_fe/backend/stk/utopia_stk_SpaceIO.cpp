#include "utopia_stk_SpaceIO.hpp"

#include "utopia_stk_MeshIO.hpp"

#include "utopia_stk_DofMap.hpp"

#include "utopia_stk_Commons.hpp"

#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace utopia {
    namespace stk {

        class SpaceIO::Impl {
        public:
            Impl(FunctionSpace &space) : space(space), io(space.mesh()) {}
            FunctionSpace &space;
            MeshIO io;
            bool is_open{false};
        };

        void SpaceIO::read(Input &in) {
            // TODO [faw 2026-07-30]: The code below was taken from
            // utopia::stk::FunctionSpace::read. The code related to
            // import_all_field_data is incomplete and needs revision.  In the
            // fluya integration the import is handled as follows:
            //
            // if (import_field_data_) {
            //     io_->import_all_field_data(true);
            //     io_->enable_interpolation_mode();
            //
            //     valid_ = io_->open_input(in);
            //     io_->load_time_step(stime_->get());
            // } else {
            //     io_->read(in);
            //     valid_ = true;
            // }
            //
            // when the todo is cleared it should be possible to just call
            // io_->read(in) from application code.
            //
            // The related code is currently commented out:
            ////////////////////////////////////////////////////////////////////
            // bool import_all_field_data = false;
            // in.get("import_all_field_data", import_all_field_data);
            //
            // if (import_all_field_data) {
            //     impl_->io.import_all_field_data(true);
            //     in.get("mesh", impl_->io);
            //
            //     if (!impl_->io.load()) {
            //         Utopia::Abort("FunctionSpace::read(in): Unable to read input!");
            //     }
            // } else {
            //     in.get("mesh", impl_->space.mesh());
            // }
            ////////////////////////////////////////////////////////////////////
            // The in.get("mesh", impl_->space.mesh()); call in the esle branch
            // does not work as it calls (for a mesh object):
            // void YAMLInput::get(const std::string &key, Configurable &val) {
            //     auto &node = impl_->root();
            //     if (node[key]) {
            //         YAMLInput branch(node[key]);
            //         val.read(branch);
            //     }
            // }
            // where val.read(branch) calls utopia::stk::Mesh::read which again
            // creates a scope local MeshIO io(*this); object.

            in.get("mesh", impl_->io);
            if (!open_input()) {
                Utopia::Abort("SpaceIO::read(in): Unable to read input!");
            }

            if (impl_->space.mesh().empty()) return;

            impl_->space.read_meta(in);
            impl_->space.register_variables();
            impl_->space.dof_map().init(impl_->space.mesh(), false);

            bool verbose = false;
            in.get("verbose", verbose);
            if (verbose) {
                std::stringstream ss;
                impl_->space.describe(ss);
                impl_->space.comm().synched_print(ss.str());
            }
        }

        bool SpaceIO::read_with_state(Input &in, Field<FunctionSpace> &field) {
            return impl_->space.read_with_state(in, field);
        }

        bool SpaceIO::write(const Field<FunctionSpace> &field) {
            if (!impl_->io.ensure_output()) {
                return false;
            }

            field.space()->backend_set_elemental_field(field);
            register_output_field(field.name());

            return impl_->io.write(1, 1);
        }

        bool SpaceIO::read_nodal(Field<FunctionSpace> &field, const bool fail_if_not_found) {
            using Bucket_t = ::stk::mesh::Bucket;
            using Communicator = Traits<FunctionSpace>::Communicator;

            auto space = field.space();
            assert(space != nullptr);

            if (!space) {
                Utopia::Abort("Space in field must be set from outside for now!");
            }

            auto &&mesh = impl_->space.mesh();
            auto &&bulk_data = mesh.bulk_data();
            auto &&meta_data = mesh.meta_data();
            auto &node_buckets = utopia::stk::universal_nodes(bulk_data);
            ::stk::mesh::FieldBase *stk_field = ::stk::mesh::get_field_by_name(field.name(), meta_data);

            if (!stk_field) {
                if (!space->comm().rank()) {
                    utopia::err() << "SpaceIO: could not find stk field with name \"" << field.name() << "\"\n";
                }

                if (fail_if_not_found) {
                    Utopia::Abort("Fatal error!");
                } else {
                    return false;
                }
            }

            int n_comp = -1;
            Vector local_vector;
            for (const auto &ib : node_buckets) {
                auto &&b = *ib;
                n_comp = ::stk::mesh::field_scalars_per_entity(*stk_field, b);

                // Resize vector
                SizeType n_nodes = utopia::stk::count_universal_nodes(bulk_data);
                SizeType nn = n_nodes * n_comp;
                local_vector.zeros(layout(Communicator::self(), nn, nn));
                local_vector.set_block_size(n_comp);
                break;
            }

            // Make sure that even if we do not have local nodes we still know how many components we have
            n_comp = space->comm().max(n_comp);

            if (!local_vector.empty()) {
                auto lv_view = local_view_device(local_vector);

                for (const auto &ib : node_buckets) {
                    auto &&b = *ib;
                    const Bucket_t::size_type length = b.size();

                    for (Bucket_t::size_type k = 0; k < length; ++k) {
                        auto node = b[k];
                        auto idx = utopia::stk::convert_entity_to_index(node);

                        if (bulk_data.in_receive_ghost(node)) {
                            idx = impl_->space.dof_map().shift_aura_idx(idx);
                        }

                        const Scalar *values = (const Scalar *)::stk::mesh::field_data(*stk_field, node);

                        for (int d = 0; d < n_comp; ++d) {
                            lv_view.set(idx * n_comp + d, values[d]);
                        }
                    }
                }
            }

            auto field_data =
                std::make_shared<Vector>(layout(space->comm(), mesh.n_local_nodes() * n_comp, mesh.n_nodes() * n_comp));
            field_data->set_block_size(n_comp);

            space->local_to_global(local_vector, *field_data, OVERWRITE_MODE);

            field.set_data(field_data);
            field.set_tensor_size(n_comp);
            field.set_space(space);

            return true;
        }

        bool SpaceIO::write(const Vector &v) { return impl_->space.write(impl_->io.output_path(), v); }

        bool SpaceIO::write(const Vector &v, const int step, const Scalar t) {
            if (!impl_->is_open) {
                open_output();
            }

            impl_->space.global_vector_to_nodal_field(v);
            return impl_->io.write(step, t);
        }

        void SpaceIO::set_output_path(const Path &path) { impl_->io.set_output_path(path); }

        void SpaceIO::set_output_mode(const std::string &output_mode) { impl_->io.set_output_mode(output_mode); }
        void SpaceIO::set_output_mode(enum OutputMode output_mode) { impl_->io.set_output_mode(output_mode); }

        bool SpaceIO::open_output() {
            if (!impl_->io.ensure_output()) {
                return false;
            }

            impl_->space.register_output_variables(impl_->io);
            impl_->is_open = true;
            return true;
        }

        SpaceIO::SpaceIO(FunctionSpace &space) : impl_(utopia::make_unique<Impl>(space)) {}
        SpaceIO::~SpaceIO() = default;

        void SpaceIO::enable_interpolation_mode() {
            assert(impl_);
            impl_->io.enable_interpolation_mode();
        }

        bool SpaceIO::load_time_step(const Scalar t) {
            assert(impl_);
            return impl_->io.load_time_step(t);
        }

        bool SpaceIO::load_last_time_step() {
            assert(impl_);
            return impl_->io.load_last_time_step();
        }

        std::vector<double> SpaceIO::get_time_steps() const { return impl_->io.get_time_steps(); }

        int SpaceIO::num_time_steps() const { return impl_->io.num_time_steps(); }
        SpaceIO::Scalar SpaceIO::max_time() const { return impl_->io.max_time(); }

        void SpaceIO::set_import_all_data(const bool val) { impl_->io.set_import_all_data(val); }

        bool SpaceIO::write(const int step, const Scalar t) {
            assert(impl_);
            return impl_->io.write(step, t);
        }

        void SpaceIO::import_all_field_data(const bool value) {
            assert(impl_);
            impl_->io.import_all_field_data(value);
        }

        void SpaceIO::set_restart_cycle_count(const int count) {
            assert(impl_);
            impl_->io.set_restart_cycle_count(count);
        }

        void SpaceIO::set_domain_decomposition_method(const std::string method) {
            assert(impl_);
            impl_->io.set_domain_decomposition_method(method);
        }

        void SpaceIO::set_input_path(const Path &path) {
            assert(impl_);
            impl_->io.set_read_path(path);
        }
        bool SpaceIO::open_input() {
            assert(impl_);
            return impl_->io.ensure_input();
        }

        bool SpaceIO::open_input(Input &in) {
            assert(impl_);
            in.get("mesh", impl_->io);
            impl_->space.read_meta(in);
            // impl_->space.register_variables();
            if (!open_input()) {
                return false;
            }

            impl_->space.dof_map().init(impl_->space.mesh());

            return true;
        }

        bool SpaceIO::open_restart(Input &in) {
            assert(impl_);
            impl_->space.read_meta(in);
            in.get("mesh", [this](Input &mesh_node) { this->impl_->io.read_minimal(mesh_node); });
            impl_->io.set_read_purpose(::stk::io::READ_RESTART);
            if (!open_input()) {
                return false;
            }

            impl_->space.dof_map().init(impl_->space.mesh());

            return true;
        }

        void SpaceIO::register_output_field(const std::string &var_name) { impl_->io.register_output_field(var_name); }

        void SpaceIO::register_output_field(const Field<FunctionSpace> &field) {
            this->update_output_field(field);
            impl_->io.register_output_field(field.name());
        }

        void SpaceIO::update_output_field(const Field<FunctionSpace> &field) {
            impl_->space.backend_set_nodal_field(field);
        }

    }  // namespace stk
}  // namespace utopia
