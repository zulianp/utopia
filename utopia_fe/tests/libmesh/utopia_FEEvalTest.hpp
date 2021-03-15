#ifndef UTOPIA_FE_EVAL_TEST_HPP
#define UTOPIA_FE_EVAL_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    class FEEvalTest final : public FETest {
    public:
        void run(Input &in) override;  // { assert(false); }

        inline static std::string command() { return "fe_test"; }
    };

}  // namespace utopia

#endif  // UTOPIA_FE_EVAL_TEST_HPP
