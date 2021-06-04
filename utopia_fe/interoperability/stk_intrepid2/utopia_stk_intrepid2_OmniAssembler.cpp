#include "utopia_stk_intrepid2_OmniAssembler.hpp"
#include "utopia_intrepid2_OmniAssembler_impl.hpp"

// Utopia::Intrepid2
#include "utopia_intrepid2_OmniAssembler.hpp"

#include "utopia_stk_intrepid2.hpp"
#include "utopia_stk_intrepid2_Transport.hpp"

namespace utopia {

    namespace intrepid2 {
        template class OmniAssembler<utopia::stk::FunctionSpace>;
    }

    class OmniAssembler<utopia::stk::FunctionSpace>::Impl {
    public:
        using Intrepid2OmniAssembler = utopia::intrepid2::OmniAssembler<utopia::stk::FunctionSpace>;

        std::unique_ptr<Intrepid2OmniAssembler> assembler;
    };

    OmniAssembler<utopia::stk::FunctionSpace>::~OmniAssembler() = default;
    OmniAssembler<utopia::stk::FunctionSpace>::OmniAssembler(const std::shared_ptr<FunctionSpace> &space)
        : impl_(utopia::make_unique<Impl>()) {
        impl_->assembler = utopia::make_unique<Impl::Intrepid2OmniAssembler>(space);
    }

    bool OmniAssembler<utopia::stk::FunctionSpace>::is_linear() const { return impl_->assembler->is_linear(); }

    void OmniAssembler<utopia::stk::FunctionSpace>::clear() { return impl_->assembler->clear(); }

    std::string OmniAssembler<utopia::stk::FunctionSpace>::name() const { return impl_->assembler->name(); }

    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(const Vector &x, Matrix &hessian, Vector &gradient) {
        return impl_->assembler->assemble(x, hessian, gradient);
    }

    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(const Vector &x, Matrix &hessian) {
        return impl_->assembler->assemble(x, hessian);
    }
    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(const Vector &x, Vector &gradient) {
        return impl_->assembler->assemble(x, gradient);
    }

    // For linear only
    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(Matrix &hessian) {
        return impl_->assembler->assemble(hessian);
    }
    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(Vector &gradient) {
        return impl_->assembler->assemble(gradient);
    }

    void OmniAssembler<utopia::stk::FunctionSpace>::read(Input &in) { impl_->assembler->read(in); }

    void OmniAssembler<utopia::stk::FunctionSpace>::set_environment(const std::shared_ptr<Environment> &env) {
        impl_->assembler->set_environment(env);
    }

    std::shared_ptr<OmniAssembler<utopia::stk::FunctionSpace>::Environment>
    OmniAssembler<utopia::stk::FunctionSpace>::environment() const {
        return impl_->assembler->environment();
    }

    void OmniAssembler<utopia::stk::FunctionSpace>::set_space(const std::shared_ptr<FunctionSpace> &space) {
        impl_->assembler->set_space(space);
    }

    std::shared_ptr<utopia::stk::FunctionSpace> OmniAssembler<utopia::stk::FunctionSpace>::space() const {
        return impl_->assembler->space();
    }

}  // namespace utopia
