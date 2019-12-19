#include "utopia.hpp"
#include "utopia_AppRunner.hpp"

int main(const int argc, char *argv[])
{
    using namespace utopia;

    Utopia::Init(argc, argv);

    {
        AppRunner app_runner;
        app_runner.run(argc, argv);
    }

    return Utopia::Finalize();
}
