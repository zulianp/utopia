
#include "utopia.hpp"
#include "utopia_Testing.hpp"

using namespace std;

int main(const int argc, char *argv[]) {
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
