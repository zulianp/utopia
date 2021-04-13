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

    }  // namespace stk
}  // namespace utopia
