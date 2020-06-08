#include "utopia.hpp"
#include "utopia_Testing.hpp"

int main(const int argc, char *argv[]) {
    using namespace utopia;

    Utopia::Init(argc, argv);
    UTOPIA_TRACE_REGION_BEGIN("main");

    {
        TestRunner app_runner;
        app_runner.run(argc, argv);
    }

    UTOPIA_TRACE_REGION_END("main");
    return Utopia::Finalize();
}
