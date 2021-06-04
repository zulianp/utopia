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

        std::shared_ptr<FunctionSpace> space;
        std::unique_ptr<FEAssembler<utopia::stk::FunctionSpace>> assembler;
    };

    OmniAssembler<utopia::stk::FunctionSpace>::~OmniAssembler() = default;
    OmniAssembler<utopia::stk::FunctionSpace>::OmniAssembler(const std::shared_ptr<FunctionSpace> &space)
        : impl_(utopia::make_unique<Impl>()) {
        impl_->space = space;
    }

    bool OmniAssembler<utopia::stk::FunctionSpace>::is_linear() const { return impl_assembler().is_linear(); }

    void OmniAssembler<utopia::stk::FunctionSpace>::clear() { return impl_assembler().clear(); }

    std::string OmniAssembler<utopia::stk::FunctionSpace>::name() const { return impl_assembler().name(); }

    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(const Vector &x, Matrix &hessian, Vector &gradient) {
        return impl_assembler().assemble(x, hessian, gradient);
    }

    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(const Vector &x, Matrix &hessian) {
        return impl_assembler().assemble(x, hessian);
    }
    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(const Vector &x, Vector &gradient) {
        return impl_assembler().assemble(x, gradient);
    }

    // For linear only
    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(Matrix &hessian) {
        return impl_assembler().assemble(hessian);
    }
    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(Vector &gradient) {
        return impl_assembler().assemble(gradient);
    }

    void OmniAssembler<utopia::stk::FunctionSpace>::read(Input &in) {
        impl_->assembler = utopia::make_unique<Impl::Intrepid2OmniAssembler>(this->space());
        impl_assembler().read(in);
    }

    void OmniAssembler<utopia::stk::FunctionSpace>::set_environment(const std::shared_ptr<Environment> &env) {
        impl_assembler().set_environment(env);
    }

    std::shared_ptr<OmniAssembler<utopia::stk::FunctionSpace>::Environment>
    OmniAssembler<utopia::stk::FunctionSpace>::environment() const {
        return impl_assembler().environment();
    }

    void OmniAssembler<utopia::stk::FunctionSpace>::set_space(const std::shared_ptr<FunctionSpace> &space) {
        impl_->space = space;
        if (impl_->assembler) {
            impl_->assembler->set_space(space);
        }
    }

    std::shared_ptr<utopia::stk::FunctionSpace> OmniAssembler<utopia::stk::FunctionSpace>::space() const {
        return impl_->space;
    }

    // Convenience function for checking not null
    OmniAssembler<utopia::stk::FunctionSpace>::ImplAssembler &
    OmniAssembler<utopia::stk::FunctionSpace>::impl_assembler() const {
        assert(impl_);
        assert(impl_->assembler);
        if (!impl_->assembler) {
            Utopia::Abort("OmniAssembler<utopia::stk::FunctionSpace> was not initialized before use!");
        }
        return *impl_->assembler;
    }

}  // namespace utopia
