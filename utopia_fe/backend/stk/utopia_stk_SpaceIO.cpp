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

        void SpaceIO::read(Input &in) { impl_->space.read(in); }

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

        bool SpaceIO::read_nodal(Field<FunctionSpace> &field) {
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

        bool SpaceIO::write(const int step, const Scalar t) {
            assert(impl_);
            return impl_->io.write(step, t);
        }

        void SpaceIO::import_all_field_data(const bool value) {
            assert(impl_);
            impl_->io.import_all_field_data(value);
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
