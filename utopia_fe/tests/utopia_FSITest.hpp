#ifndef UTOPIA_FSI_TEST_HPP
#define UTOPIA_FSI_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    class FSITest final : public FETest {
    public:
        void run(Input &in) override;  // { assert(false && "IMPLEMENT ME"); }

        inline static std::string command() { return "fsi"; }
    };

}  // namespace utopia

#endif  // UTOPIA_FSI_TEST_HPP
