#include "utopia_Main.hpp"

namespace utopia {

    int Main(const int argc, char *argv[]) {
        using namespace utopia;

        Utopia::Init(argc, argv);
        UTOPIA_TRACE_REGION_BEGIN("main");

        {
            AppRunner app_runner;
            app_runner.run(argc, argv);
        }

        UTOPIA_TRACE_REGION_END("main");
        return Utopia::Finalize();
    }
}  // namespace utopia
