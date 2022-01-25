#ifndef UTOPIA_APP_REGISTRY_H
#define UTOPIA_APP_REGISTRY_H

#include <functional>
#include <iostream>
#include <map>
#include <string>

#include "utopia_ActionRegistry.hpp"

namespace utopia {

    class AppRegistry {
    public:
        using Count = long;

        using RunApp = void (*)();

        static char add_app(const std::string &app_name, RunApp run_app);
        static char add_app(const std::string &app_name, const std::string &backend, RunApp run_app);

        static AppRegistry &instance();
        int run(const std::string &app_name);
        int run(const std::string &app_name, const std::string &backend);
        int run_all();
        void describe(std::ostream &os = std::cout) const;
        bool verbose() const;
        void verbose(const bool val);

    private:
        AppRegistry();
        ~AppRegistry();
        ActionRegistry apps_;

        static std::string concat_app_backend(const std::string &app_name, const std::string &backend);
    };

}  // namespace utopia

#endif  // UTOPIA_APP_REGISTRY_H
