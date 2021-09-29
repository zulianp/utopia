#include "utopia_mars_FEAssembler.hpp"
#include "utopia_mars_Factory.hpp"

#include "utopia_mars_LaplaceOperator.hpp"
#include "utopia_mars_LinearElasticity.hpp"
#include "utopia_mars_Mass.hpp"
#include "utopia_mars_WeakLinearThermoElasticity.hpp"

namespace utopia {
    namespace mars {

        template <class... Args>
        class ConcreteFactory final : public Factory {
        public:
            inline std::unique_ptr<FEAssembler> new_assembler(Input &in) const override {
                std::string type;
                in.get("type", type);
                return new_assembler(type);
            }

            inline std::unique_ptr<FEAssembler> new_assembler(const std::string &type) const {
                if (type == "LaplaceOperator") {
                    return utopia::make_unique<LaplaceOperator<Args...>>();
                } else if (type == "LinearElasticity") {
                    return utopia::make_unique<LinearElasticity<Args...>>();
                } else if (type == "Mass") {
                    return utopia::make_unique<Mass<Args...>>();
                } else if (type == "WeakLinearThermoElasticity") {
                    return utopia::make_unique<WeakLinearThermoElasticity<Args...>>();
                } else {
                    return nullptr;
                }
            }

            inline static ConcreteFactory &instance() {
                static ConcreteFactory instance_;
                return instance_;
            }

        private:
            ConcreteFactory() = default;

            // Add more mars object
        };

    }  // namespace mars
}  // namespace utopia
