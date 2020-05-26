#ifndef UTOPIA_SEMIGEOMETRIC_MULTIGRID_TEST_HPP
#define UTOPIA_SEMIGEOMETRIC_MULTIGRID_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    class SMGTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "smg"; }
    };

}  // namespace utopia

#endif  // UTOPIA_SEMIGEOMETRIC_MULTIGRID_TEST_HPP
