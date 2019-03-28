#ifndef UTOPIA_NON_LINEAR_ELASTICITY_TEST_HPP
#define UTOPIA_NON_LINEAR_ELASTICITY_TEST_HPP

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

    class NonLinearElasticityTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "nlelast";
        }
    };

}

#endif //UTOPIA_NON_LINEAR_ELASTICITY_TEST_HPP
