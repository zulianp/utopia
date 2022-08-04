#include "utopia_stk_intrepid2_OmniAssembler.hpp"
#include "utopia_kokkos_OmniAssembler_impl.hpp"

// Utopia::Intrepid2
#include "utopia_intrepid2_OmniAssembler.hpp"

#include "utopia_stk_intrepid2.hpp"
#include "utopia_stk_intrepid2_L2Projection.hpp"
#include "utopia_stk_intrepid2_Transport.hpp"

namespace utopia {

    namespace kokkos {
        template class OmniAssembler<utopia::stk::FunctionSpace,
                                     utopia::intrepid2::FE<utopia::stk::FunctionSpace::Scalar>>;
    }

    class OmniAssembler<utopia::stk::FunctionSpace>::Impl {
    public:
        using Intrepid2OmniAssembler = utopia::intrepid2::OmniAssembler<utopia::stk::FunctionSpace>;
        using Intrepid2FE = utopia::intrepid2::FE<Scalar>;

        std::shared_ptr<Environment> environment;
        std::shared_ptr<FunctionSpace> space;
        std::unique_ptr<Intrepid2OmniAssembler> assembler;
        std::shared_ptr<Intrepid2FE> fe;

        Intrepid2OmniAssembler &get_assembler() {
            assert(assembler);
            if (!assembler) {
                Utopia::Abort("OmniAssembler<utopia::stk::FunctionSpace> was not initialized before use!");
            }
            return *assembler;
        }
    };

    OmniAssembler<utopia::stk::FunctionSpace>::~OmniAssembler() = default;
    OmniAssembler<utopia::stk::FunctionSpace>::OmniAssembler(const std::shared_ptr<FunctionSpace> &space)
        : impl_(utopia::make_unique<Impl>()) {
        impl_->space = space;
    }

    void OmniAssembler<utopia::stk::FunctionSpace>::set_time(const std::shared_ptr<SimulationTime> &time) {
        assert(impl_->assembler);
        if (impl_->assembler) {
            impl_->assembler->set_time(time);
        }
    }

    bool OmniAssembler<utopia::stk::FunctionSpace>::is_linear() const { return impl_->get_assembler().is_linear(); }
    bool OmniAssembler<utopia::stk::FunctionSpace>::is_operator() const { return impl_->get_assembler().is_operator(); }

    void OmniAssembler<utopia::stk::FunctionSpace>::clear() { return impl_->get_assembler().clear(); }

    std::string OmniAssembler<utopia::stk::FunctionSpace>::name() const { return impl_->get_assembler().name(); }

    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(const Vector &x, Matrix &hessian, Vector &gradient) {
        return impl_->get_assembler().assemble(x, hessian, gradient);
    }

    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(const Vector &x, Matrix &hessian) {
        return impl_->get_assembler().assemble(x, hessian);
    }
    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(const Vector &x, Vector &gradient) {
        return impl_->get_assembler().assemble(x, gradient);
    }

    bool OmniAssembler<utopia::stk::FunctionSpace>::apply(const Vector &x, Vector &hessian_times_x) {
        return impl_->get_assembler().apply(x, hessian_times_x);
    }

    // For linear only
    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(Matrix &hessian) {
        return impl_->get_assembler().assemble(hessian);
    }
    bool OmniAssembler<utopia::stk::FunctionSpace>::assemble(Vector &gradient) {
        return impl_->get_assembler().assemble(gradient);
    }

    void OmniAssembler<utopia::stk::FunctionSpace>::read(Input &in) {
        impl_->assembler = nullptr;

        if (!impl_->fe) {
            int quadrature_order = 2;
            in.get("quadrature_order", quadrature_order);
            impl_->fe = std::make_shared<Impl::Intrepid2FE>();
            create_fe(*impl_->space, *impl_->fe, quadrature_order);
        }

        if (!impl_->assembler) {
            impl_->assembler = utopia::make_unique<Impl::Intrepid2OmniAssembler>(this->space());
            impl_->assembler->set_domain_fe(impl_->fe);
            impl_->assembler->set_environment(impl_->environment);
        }

        in.get("material", [this](Input &node) {
            std::string type;
            node.get("type", type);

            if (!type.empty()) {
                if (type == "Transport") {
                    auto ass = utopia::make_unique<stk::Transport>(impl_->fe);
                    // ass->set_fe(impl_->fe);
                    ass->set_space(this->space());
                    ass->set_environment(impl_->environment);
                    ass->read(node);

                    // We add it from the outside
                    impl_->get_assembler().fail_if_unregistered(false);
                    impl_->get_assembler().add_domain_assembler(std::move(ass));
                } else if (type == "L2Projection") {
                    auto ass = utopia::make_unique<stk::L2Projection>(impl_->fe);
                    // ass->set_fe(impl_->fe);
                    ass->set_space(this->space());
                    ass->set_environment(impl_->environment);
                    ass->read(node);

                    // We add it from the outside
                    impl_->get_assembler().fail_if_unregistered(false);
                    impl_->get_assembler().add_domain_assembler(std::move(ass));
                }
            }
        });

        impl_->get_assembler().read(in);
    }

    void OmniAssembler<utopia::stk::FunctionSpace>::set_environment(const std::shared_ptr<Environment> &env) {
        impl_->environment = env;

        if (impl_->assembler) {
            impl_->assembler->set_environment(env);
        }
    }

    std::shared_ptr<OmniAssembler<utopia::stk::FunctionSpace>::Environment>
    OmniAssembler<utopia::stk::FunctionSpace>::environment() const {
        return impl_->get_assembler().environment();
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

}  // namespace utopia
