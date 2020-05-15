#ifndef UTOPIA_FE_REFACTORED_CONTACT_TEST_HPP
#define UTOPIA_FE_REFACTORED_CONTACT_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    class RefactoredContactTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "rcontact"; }
    };

}  // namespace utopia

#endif  // UTOPIA_FE_REFACTORED_CONTACT_TEST_HPP
