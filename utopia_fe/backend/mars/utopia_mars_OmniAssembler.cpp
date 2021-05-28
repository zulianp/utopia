#include "utopia_mars_OmniAssembler.hpp"

namespace utopia {
    namespace mars {

        class OmniAssembler::Impl {
        public:
        };

        OmniAssembler::OmniAssembler(const std::shared_ptr<mars::FunctionSpace> &space) {}

        OmniAssembler::~OmniAssembler() {}

        bool OmniAssembler::assemble(const Vector &x, Matrix &jacobian, Vector &fun) {}

        bool OmniAssembler::assemble(const Vector &x, Matrix &jacobian) {}

        bool OmniAssembler::assemble(const Vector &x, Vector &fun) {}

        bool OmniAssembler::assemble(Matrix &jacobian) {}

        void OmniAssembler::read(Input &in) {}

        void OmniAssembler::set_environment(const std::shared_ptr<Environment<mars::FunctionSpace>> &env) {}

        bool OmniAssembler::is_linear() const {}

    }  // namespace mars
}  // namespace utopia
