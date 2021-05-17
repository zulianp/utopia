#ifndef UTOPIA_MECH_TEST_HPP
#define UTOPIA_MECH_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    class MechTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "mech"; }
    };

}  // namespace utopia

#endif  // UTOPIA_MECH_TEST_HPP
