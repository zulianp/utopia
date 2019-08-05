#ifndef UTOPIA_FE_ENERGY_ASSEMBLY_TEST_HPP
#define UTOPIA_FE_ENERGY_ASSEMBLY_TEST_HPP

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

    class EnergyTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "energy";
        }
    };

}



#endif //UTOPIA_FE_ENERGY_ASSEMBLY_TEST_HPP
