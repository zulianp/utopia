#ifndef UTOPIA_REFINE_APP_HPP
#define UTOPIA_REFINE_APP_HPP

#include <string>
#include "utopia_FEApp.hpp"
#include "utopia_UIForcingFunction.hpp"

namespace utopia {
    class RefineApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "-refine"; }
    };

}  // namespace utopia

#endif  // UTOPIA_POISSON_APP_HPP
