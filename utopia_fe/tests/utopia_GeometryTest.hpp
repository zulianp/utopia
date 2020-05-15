#ifndef UTOPIA_GEOMETRY_TEST_HPP
#define UTOPIA_GEOMETRY_TEST_HPP

#include <string>
#include "utopia_FETest.hpp"
#include "utopia_fe_base.hpp"

namespace utopia {

    class GeometryTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "geo"; }
    };

}  // namespace utopia

#endif  // UTOPIA_GEOMETRY_TEST_HPP
