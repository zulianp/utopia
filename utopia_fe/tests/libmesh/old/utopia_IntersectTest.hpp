#ifndef UTOPIA_INTERSECT_TEST_HPP
#define UTOPIA_INTERSECT_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    class IntersectTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "isect"; }
    };

}  // namespace utopia

#endif  // UTOPIA_INTERSECT_TEST_HPP
