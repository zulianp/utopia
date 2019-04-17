#ifndef UTOPIA_SEMIGEOMETRIC_MULTIGRID_TEST_HPP
#define UTOPIA_SEMIGEOMETRIC_MULTIGRID_TEST_HPP


#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

    class SMGTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "smg";
        }
    };

}


#endif //UTOPIA_SEMIGEOMETRIC_MULTIGRID_TEST_HPP
