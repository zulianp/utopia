#include "utopia.hpp"
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
        utopia::out() << "Hello this is UTOPIA_GIT_VERSION: " << UTOPIA_GIT_VERSION << std::endl;
    }

    // Finalizes library and dependencies runtimes
    return Utopia::Finalize();
}