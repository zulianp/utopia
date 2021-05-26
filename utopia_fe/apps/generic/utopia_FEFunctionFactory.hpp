#ifndef UTOPIA_FE_FUNCTION_FACTORY_HPP
#define UTOPIA_FE_FUNCTION_FACTORY_HPP

#include "utopia_FEModelFunction.hpp"
#include "utopia_ImplicitEulerIntegrator.hpp"
#include "utopia_NewmarkIntegrator.hpp"

namespace utopia {

    template <class FunctionSpace>
    class CoupledFEFunction;

    template <class FunctionSpace>
    class FEFunctionFactory {
    public:
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;
        using ImplicitEulerIntegrator_t = utopia::ImplicitEulerIntegrator<FunctionSpace>;
        using NewmarkIntegrator_t = utopia::NewmarkIntegrator<FunctionSpace>;
        using TimeDependentFunction_t = utopia::TimeDependentFunction<FunctionSpace>;

        static std::unique_ptr<FEFunctionInterface<FunctionSpace>> make(const std::shared_ptr<FunctionSpace> &space,
                                                                        Input &in) {
            std::string integrator;
            in.get("integrator", integrator);
            auto problem = make(space, integrator);
            return problem;
        }

        static std::unique_ptr<FEFunctionInterface<FunctionSpace>> make(const std::shared_ptr<FunctionSpace> &space,
                                                                        const std::string &integrator) {
            if (integrator == "Newmark") {
                return utopia::make_unique<NewmarkIntegrator_t>(space);
            } else if (integrator == "ImplicitEuler") {
                return utopia::make_unique<ImplicitEulerIntegrator_t>(space);
            } else {
                return std::make_shared<FEModelFunction_t>(space);
            }
        }

        static inline FEFunctionFactory &instance() {
            static FEFunctionFactory instance_;
            return instance_;
        }

    private:
        FEFunctionFactory() = default;
    };

}  // namespace utopia

#endif  // UTOPIA_FE_FUNCTION_FACTORY_HPP