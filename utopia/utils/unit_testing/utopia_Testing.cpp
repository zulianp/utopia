#include "utopia_Testing.hpp"
#include "utopia_Tracer.hpp"

namespace utopia {
    int Test(const int argc, char *argv[]) {
        using namespace utopia;

        Utopia::Init(argc, argv);
        UTOPIA_TRACE_REGION_BEGIN("main");

        {
            TestRunner test_runner;
            test_runner.run(argc, argv);
        }

        UTOPIA_TRACE_REGION_END("main");
        return Utopia::Finalize();
    }

}  // namespace utopia