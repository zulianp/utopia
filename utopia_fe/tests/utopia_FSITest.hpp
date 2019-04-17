#ifndef UTOPIA_FSI_TEST_HPP
#define UTOPIA_FSI_TEST_HPP

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

    class FSITest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "fsi";
        }
    };

}

#endif //UTOPIA_FSI_TEST_HPP
