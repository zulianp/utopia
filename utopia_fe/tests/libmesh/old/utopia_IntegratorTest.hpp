#ifndef UTOPIA_INTEGRATOR_TEST_HPP
#define UTOPIA_INTEGRATOR_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    class IntegratorTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "integrator"; }
    };

}  // namespace utopia

#endif  // UTOPIA_INTEGRATOR_TEST_HPP
