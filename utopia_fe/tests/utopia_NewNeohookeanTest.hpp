#ifndef UTOPIA_NEW_NEOHOKEAN_TEST_HPP
#define UTOPIA_NEW_NEOHOKEAN_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"

namespace utopia {

    class NewNeohookeanTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "neohook"; }
    };

}  // namespace utopia

#endif  // UTOPIA_NEW_NEOHOKEAN_TEST_HPP
