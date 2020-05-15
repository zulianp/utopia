#ifndef UTOPIA_MSH_READER_TEST_HPP
#define UTOPIA_MSH_READER_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    class MSHReaderTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "msh"; }
    };

}  // namespace utopia

#endif  // UTOPIA_MSH_READER_TEST_HPP
