#include "utopia.hpp"
#include "utopia_Testing.hpp"
#include "utopia_fe_Instance.hpp"

int main(int argc, char *argv[]) {
    using namespace utopia;

    UtopiaFE::Init(argc, argv);
    UTOPIA_TRACE_REGION_BEGIN("main");

    {
        TestRunner app_runner;
        app_runner.run(argc, argv);
    }

    UTOPIA_TRACE_REGION_END("main");
    return UtopiaFE::Finalize();
}
