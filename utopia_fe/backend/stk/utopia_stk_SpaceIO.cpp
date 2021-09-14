#include "utopia_stk_SpaceIO.hpp"

#include "utopia_stk_MeshIO.hpp"

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

            return true;
        }

        void SpaceIO::register_output_field(const std::string &var_name) { impl_->io.register_output_field(var_name); }

    }  // namespace stk
}  // namespace utopia
