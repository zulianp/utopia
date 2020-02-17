#ifndef UTOPIA_DUAL_BASIS_TEST_HPP
#define UTOPIA_DUAL_BASIS_TEST_HPP

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

    class DualBasisTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "dualbasis";
        }
    };

}

#endif //UTOPIA_DUAL_BASIS_TEST_HPP
