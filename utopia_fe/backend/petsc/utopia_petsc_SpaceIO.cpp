#include "utopia_petsc_SpaceIO.hpp"
#include "utopia_petsc_IO.hpp"

#include "utopia_petsc_DMDABase.hpp"

namespace utopia {
    namespace petsc {

        class SpaceIO::Impl {
        public:
            Impl(FunctionSpace &space) : space(space) /*, io(space.mesh())*/ {}
            FunctionSpace &space;
            // PetscIO io;
            Path output_path{"out.vtr"};
        };

        void SpaceIO::read(Input &in) { impl_->space.read(in); }

        bool SpaceIO::read_with_state(Input &in, Field<FunctionSpace> &field) {
            // return impl_->space.read_with_state(in, field);
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
            return false;
        }

        bool SpaceIO::write(const Vector &v) { return impl_->space.write(impl_->output_path, v); }

        bool SpaceIO::write(const Vector &v, const int step, const Scalar t) {
            Path path = impl_->output_path.without_extension() + "_" + std::to_string(step) + "." + "vtr";

            if (!impl_->space.write(path, v)) {
                utopia::out().stream() << "Failed to write: " << path << "\n";
                return false;
            } else {
                utopia::out().stream() << "Wrote: " << path << "\n";
                return true;
            }
        }

        void SpaceIO::set_output_path(const Path &path) { impl_->output_path = path; }

        bool SpaceIO::open_output() { return true; }

        SpaceIO::SpaceIO(FunctionSpace &space) : impl_(utopia::make_unique<Impl>(space)) {}
        SpaceIO::~SpaceIO() = default;

        void SpaceIO::enable_interpolation_mode() {
            // assert(impl_);
            // impl_->io.enable_interpolation_mode();

            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        bool SpaceIO::load_time_step(const Scalar t) {
            // assert(impl_);
            // return impl_->io.load_time_step(t);
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
            return false;
        }

        bool SpaceIO::write(const int step, const Scalar t) {
            // assert(impl_);
            // return impl_->io.write(step, t);
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
            return false;
        }

        void SpaceIO::import_all_field_data(const bool value) {
            // assert(impl_);
            // impl_->io.import_all_field_data(value);

            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

        void SpaceIO::set_input_path(const Path &path) {
            // assert(impl_);
            // impl_->io.set_read_path(path);

            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }
        bool SpaceIO::open_input() {
            // assert(impl_);
            // return impl_->io.ensure_input();
            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
            return false;
        }

        bool SpaceIO::open_input(Input &in) {
            // assert(impl_);
            // in.get("mesh", impl_->io);
            // impl_->space.read_meta(in);
            // // impl_->space.register_variables();
            // if (!open_input()) {
            //     return false;
            // }

            // return true;

            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
            return false;
        }

        void SpaceIO::register_output_field(const std::string &var_name) {
            // impl_->io.register_output_field(var_name);

            assert(false);
            Utopia::Abort("IMPLEMENT ME!");
        }

    }  // namespace petsc
}  // namespace utopia
