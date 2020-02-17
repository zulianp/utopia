#include "utopia.hpp"
#include "utopia_benchmarks.hpp"

using namespace std;

int main(const int argc, char *argv[])
{
    using namespace utopia;

    Utopia::Init(argc, argv);

    {
        run_benchmarks();
    }

    return Utopia::Finalize();
}
