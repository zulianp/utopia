#ifndef UTOPIA_DUAL_BASIS_TEST_HPP
#define UTOPIA_DUAL_BASIS_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    class DualBasisTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "dualbasis"; }
    };

}  // namespace utopia

#endif  // UTOPIA_DUAL_BASIS_TEST_HPP
