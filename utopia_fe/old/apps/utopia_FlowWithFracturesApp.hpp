#ifndef UTOPIA_FLOW_WITH_FRACTURES_APP_HPP
#define UTOPIA_FLOW_WITH_FRACTURES_APP_HPP

#include <string>
#include "utopia_FEApp.hpp"
#include "utopia_UIForcingFunction.hpp"

namespace utopia {
    class FlowWithFracturesApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "-flow_wf"; }
    };
}  // namespace utopia

#endif  // UTOPIA_FLOW_WITH_FRACTURES_APP_HPP
