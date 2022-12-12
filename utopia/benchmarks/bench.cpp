#include "utopia.hpp"
#include "utopia_benchmarks.hpp"

using namespace std;

int main(const int argc, char *argv[]) {
    using namespace utopia;

    Utopia::Init(argc, argv);
    UTOPIA_TRACE_REGION_BEGIN("main");

    run_benchmarks();

    UTOPIA_TRACE_REGION_BEGIN("end");
    return Utopia::Finalize();
}
