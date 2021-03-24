#ifndef UTOPIA_ASSEMBLY_TEST_HPP
#define UTOPIA_ASSEMBLY_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"

namespace utopia {

    class AssemblyTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "asm"; }
    };

}  // namespace utopia

#endif  // UTOPIA_ASSEMBLY_TEST_HPP
