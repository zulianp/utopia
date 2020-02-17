#ifndef UTOPIA_FE_TENSOR_TEST_HPP
#define UTOPIA_FE_TENSOR_TEST_HPP

#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

    class FETensorTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "tensor";
        }
    };

}

#endif //UTOPIA_FE_TENSOR_TEST_HPP
