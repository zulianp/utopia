#ifndef UTOPIA_NON_LINEAR_ELASTICITY_TEST_HPP
#define UTOPIA_NON_LINEAR_ELASTICITY_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    class NonLinearElasticityTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "nlelast"; }
    };

}  // namespace utopia

#endif  // UTOPIA_NON_LINEAR_ELASTICITY_TEST_HPP
