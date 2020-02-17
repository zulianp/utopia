#ifndef UTOPIA_INTEGRATOR_TEST_HPP
#define UTOPIA_INTEGRATOR_TEST_HPP

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

    class IntegratorTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "integrator";
        }
    };

}

#endif //UTOPIA_INTEGRATOR_TEST_HPP
