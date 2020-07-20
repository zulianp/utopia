#include "utopia.hpp"
#include "utopia_Reporter.hpp"
#include "utopia_Version.hpp"

#include <iostream>

// Run it with `./examples/example_hello_world
int main(const int argc, char *argv[]) {
    using namespace utopia;

    // Initializes library and dependencies runtimes
    Utopia::Init(argc, argv);

    {
        // Scope for utopia programs if written in the main function
        // Your utopia code here...

        // Print git version to terminal
        out() << "Hello this is UTOPIA_GIT_VERSION: " << UTOPIA_GIT_VERSION << '\n';

        out() << "Use utopia::out() instead of std::cout!\n";
        err() << "Use utopia::err() instead of std::cerr!\n";
        dev() << "Use utopia::dev() if your output is meant for developers!\n";
    }

    // Finalizes library and dependencies runtimes
    return Utopia::Finalize();
}