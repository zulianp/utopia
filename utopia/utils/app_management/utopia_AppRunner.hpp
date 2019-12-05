#ifndef UTOPIA_APP_RUNNER_HPP
#define UTOPIA_APP_RUNNER_HPP

#include <string>
#include <vector>
#include <iostream>
#include "utopia_AppMacros.hpp"

namespace utopia {

    class AppRunner {
    public:
        AppRunner();
        virtual ~AppRunner();
        int run(int argc, char **argv) const;
        int run(const std::vector<std::string> &apps) const;
        void describe(std::ostream &os = std::cout) const;
        void verbose(const bool val);
    };
}

#endif // UTOPIA_APP_RUNNER_HPP
