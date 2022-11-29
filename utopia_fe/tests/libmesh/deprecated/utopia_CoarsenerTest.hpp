#ifndef UTOPIA_COARSENER_TEST_HPP
#define UTOPIA_COARSENER_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    class CoarsenerTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "coarsener_test"; }
    };

}  // namespace utopia

#endif  // UTOPIA_COARSENER_TEST_HPP
