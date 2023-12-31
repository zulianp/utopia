#ifndef UTOPIA_FE_FUNCTION_FACTORY_HPP
#define UTOPIA_FE_FUNCTION_FACTORY_HPP

#include "utopia_fe_config.hpp"

#include "utopia_FEModelFunction.hpp"
#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_IncrementalLoading.hpp"
#include "utopia_NewmarkIntegrator.hpp"
#include "utopia_VelocityImplicitEulerIntegrator.hpp"
#include "utopia_VelocityNewmarkIntegrator.hpp"

namespace utopia {

    template <class FunctionSpace>
    class CoupledFEFunction;

    template <class FunctionSpace>
    class FEFunctionFactoryBase {
    public:
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        using ImplicitEulerIntegrator_t = utopia::ImplicitEulerIntegrator<FunctionSpace>;
        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;

        using VelocityImplicitEulerIntegrator_t = utopia::VelocityImplicitEulerIntegrator<FunctionSpace>;
        using VelocityNewmarkIntegrator_t = utopia::VelocityNewmarkIntegrator<FunctionSpace>;

        using TimeDependentFunction_t = utopia::TimeDependentFunction<FunctionSpace>;
        using IncrementalLoading_t = utopia::IncrementalLoading<FunctionSpace>;

        static std::unique_ptr<FEFunctionInterface<FunctionSpace>> make(const std::shared_ptr<FunctionSpace> &space,
                                                                        Input &in) {
            std::string integrator;
            in.get("integrator", integrator);
            auto problem = make(space, integrator);
            return problem;
        }

        static std::unique_ptr<FEFunctionInterface<FunctionSpace>> make(const std::shared_ptr<FunctionSpace> &space,
                                                                        const std::string &integrator,
                                                                        const bool fail_if_unregistered = true) {
            return make_time_integrator(
                utopia::make_unique<FEModelFunction_t>(space), integrator, fail_if_unregistered);
        }

        static std::unique_ptr<FEFunctionInterface<FunctionSpace>> make_time_integrator(
            const std::shared_ptr<FEFunctionInterface<FunctionSpace>> &fun,
            const std::string &integrator,
            const bool fail_if_unregistered = true) {
            if (integrator == "Newmark") {
                return utopia::make_unique<NewmarkIntegrator_t>(std::move(fun));
            } else if (integrator == "ImplicitEuler") {
                return utopia::make_unique<ImplicitEulerIntegrator_t>(std::move(fun));
            } else if (integrator == "VelocityNewmark") {
                return utopia::make_unique<VelocityNewmarkIntegrator_t>(std::move(fun));
            } else if (integrator == "VelocityImplicitEuler") {
                return utopia::make_unique<VelocityImplicitEulerIntegrator_t>(std::move(fun));
            } else if (integrator == "IncrementalLoading") {
                return utopia::make_unique<IncrementalLoading_t>(std::move(fun));
            } else if (fail_if_unregistered) {
                Utopia::Abort("Specified time integrator does not exists!");
            }

            return nullptr;
        }

        static inline FEFunctionFactoryBase &instance() {
            static FEFunctionFactoryBase instance_;
            return instance_;
        }

        virtual ~FEFunctionFactoryBase() = default;

    private:
        FEFunctionFactoryBase() = default;
    };

    template <class FunctionSpace>
    class FEFunctionFactory : public FEFunctionFactoryBase<FunctionSpace> {};

}  // namespace utopia

#endif  // UTOPIA_FE_FUNCTION_FACTORY_HPP