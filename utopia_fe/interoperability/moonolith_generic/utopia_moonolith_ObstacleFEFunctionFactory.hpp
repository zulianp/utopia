#ifndef UTOPIA_MOONOLITH_OBSTACLE_FE_FUNCTION_FACTORY_HPP
#define UTOPIA_MOONOLITH_OBSTACLE_FE_FUNCTION_FACTORY_HPP

#include "utopia_FEFunctionFactory.hpp"

#include "utopia_ObstacleNewmark.hpp"
#include "utopia_ObstacleStabilizedNewmark.hpp"
#include "utopia_ObstacleVelocityNewmark.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ObstacleFEFunctionFactory {
    public:
        using FEModelFunction_t = utopia::FEModelFunction<FunctionSpace>;

        using ObstacleNewmark_t = utopia::ObstacleNewmark<FunctionSpace>;
        using ObstacleVelocityNewmark_t = utopia::ObstacleVelocityNewmark<FunctionSpace>;
        using ObstacleStabilizedNewmark_t = utopia::ObstacleStabilizedNewmark<FunctionSpace>;

        static std::unique_ptr<FEFunctionInterface<FunctionSpace>> make(const std::shared_ptr<FunctionSpace> &space,
                                                                        Input &in,
                                                                        const bool fail_if_unregistered = true) {
            std::string integrator;
            in.get("integrator", integrator);
            auto problem = make(space, integrator, fail_if_unregistered);
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
            auto &&integr = FEFunctionFactoryBase<FunctionSpace>::make_time_integrator(fun, integrator, false);

            if (integr) {
                return std::move(integr);
            } else if (integrator == "ObstacleVelocityNewmark") {
                return utopia::make_unique<ObstacleVelocityNewmark_t>(std::move(fun));
            } else if (integrator == "ObstacleNewmark") {
                return utopia::make_unique<ObstacleNewmark_t>(std::move(fun));
            } else if (integrator == "ObstacleStabilizedNewmark") {
                return utopia::make_unique<ObstacleStabilizedNewmark_t>(std::move(fun));
            } else if (fail_if_unregistered) {
                Utopia::Abort("Specified time integrator does not exists!");
            }

            return nullptr;
        }

        static inline ObstacleFEFunctionFactory &instance() {
            static ObstacleFEFunctionFactory instance_;
            return instance_;
        }

        virtual ~ObstacleFEFunctionFactory() = default;

    private:
        ObstacleFEFunctionFactory() = default;
    };

}  // namespace utopia

#endif  // UTOPIA_MOONOLITH_OBSTACLE_FE_FUNCTION_FACTORY_HPP