#ifndef UTOPIA_COARSENER_TEST_HPP
#define UTOPIA_COARSENER_TEST_HPP

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

    class CoarsenerTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "coarsener_test";
        }
    };

}


#endif //UTOPIA_COARSENER_TEST_HPP
