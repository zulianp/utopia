#include "utopia_stk_SpaceIO.hpp"

#include "utopia_stk_MeshIO.hpp"

namespace utopia {
    namespace stk {

        class SpaceIO::Impl {
        public:
            Impl(FunctionSpace &space) : space(space), io(space.mesh()) {}
            FunctionSpace &space;

            // Path read_path;
            Path output_path{"out.e"};
            MeshIO io;
        };

        void SpaceIO::read(Input &in) { impl_->space.read(in); }

        bool SpaceIO::read_with_state(Input &in, Field<FunctionSpace> &field) {
            return impl_->space.read_with_state(in, field);
        }

        // bool SpaceIO::load() {}

        bool SpaceIO::write(const Vector &v) { return impl_->space.write(impl_->output_path, v); }

        bool SpaceIO::write(const Vector &v, const int step, const Scalar t) {
            assert(false && "IMPLEMENT ME");
            return false;

            // try {
            //     auto &sol = *impl_->space.raw_type_system().solution;
            //     utopia::convert(v, sol);
            //     sol.close();
            //     impl_->io.write_timestep(impl_->output_path.to_string(), impl_->space.raw_type(), step, t);
            //     return true;
            // } catch (const std::exception &ex) {
            //     utopia::err() << "[Error] SpaceIO::write: " << ex.what() << "\n";
            //     return false;
            // }
        }

        // bool SpaceIO::read(Vector &v, const int step, const Scalar t) {}

        void SpaceIO::set_output_path(const Path &path) { impl_->output_path = path; }

        // void SpaceIO::set_read_path(const Path &path) { impl_->read_path = path; }

        SpaceIO::SpaceIO(FunctionSpace &space) : impl_(utopia::make_unique<Impl>(space)) {}
        SpaceIO::~SpaceIO() = default;

    }  // namespace stk
}  // namespace utopia