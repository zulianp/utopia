#ifndef UTOPIA_MASS_APP_HPP
#define UTOPIA_MASS_APP_HPP

#include <string>
#include "utopia_FEApp.hpp"
#include "utopia_UIForcingFunction.hpp"

namespace utopia {
    class MassApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "-mass"; }

    private:
#ifdef UTOPIA_WITH_TINY_EXPR
        std::shared_ptr<SymbolicFunction> fun;
#else
        std::shared_ptr<ConstantCoefficient<double, 0> > fun;
#endif  // UTOPIA_WITH_TINY_EXPR

        bool fun_is_constant;
    };

}  // namespace utopia

#endif  // UTOPIA_MASS_APP_HPP
