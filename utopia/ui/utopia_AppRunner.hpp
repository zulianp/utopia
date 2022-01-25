#ifndef UTOPIA_APP_RUNNER_HPP
#define UTOPIA_APP_RUNNER_HPP

#include <iostream>
#include <string>
#include <vector>
#include "utopia_AppMacros.hpp"
#include "utopia_AppWithInputRegistry.hpp"
#include "utopia_Input.hpp"

namespace utopia {

    class AppRunner {
    public:
        AppRunner();
        virtual ~AppRunner();
        int run(int argc, char **argv);
        int run(const std::vector<std::string> &apps, const std::string &backend = "") const;
        inline int run(const std::string &name, Input &in) const { return run(name, "", in); }
        int run(const std::string &name, const std::string &backend, Input &in) const;
        void describe(std::ostream &os = std::cout) const;
        void verbose(const bool val);

        inline static char register_app(const std::string &name, AppRegistry::RunApp app) {
            return utopia::AppRegistry::add_app(name, app);
        }

        inline static char register_app(const std::string &name, AppWithInputRegistry::ExecuteAction app) {
            return utopia::AppWithInputRegistry::instance().add_action(name, app);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_APP_RUNNER_HPP
