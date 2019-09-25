#include "utopia_TestRunner.hpp"
#include "utopia_TestRegistry.hpp"
#include "utopia_Base.hpp"


#include <iostream>

namespace utopia {

    TestRunner::TestRunner()
    {}

    TestRunner::~TestRunner()
    {}

    int TestRunner::run(int argc, char **argv) const {
        UTOPIA_UNUSED(argc);
        UTOPIA_UNUSED(argv);
        TestRegistry::instance().describe();
        return TestRegistry::instance().run_all();
    }

}

