#ifndef UTOPIA_BOUNDARY_INTEGRAL_TEST_HPP
#define UTOPIA_BOUNDARY_INTEGRAL_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"

namespace utopia {

    class BoundaryIntegralTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "bit"; }
    };

}  // namespace utopia

#endif  // UTOPIA_BOUNDARY_INTEGRAL_TEST_HPP
