#ifndef UTOPIA_APP_REGISTRY_H
#define UTOPIA_APP_REGISTRY_H

#include <map>
#include <string>
#include <functional>
#include <iostream>

#include "utopia_ActionRegistry.hpp"

namespace utopia {

    class AppRegistry {
    public:
        using Count = long;

        typedef void (*RunApp)();

        static char add_app(const std::string &app_name, RunApp run_app);

        static AppRegistry &instance();
        int run(const std::string &app_name);
        int run_all();
        void describe(std::ostream &os = std::cout) const;
        bool verbose() const;
        void verbose(const bool val);

    private:
        AppRegistry();
        ~AppRegistry();
        ActionRegistry apps_;
    };

}

#endif //UTOPIA_APP_REGISTRY_H

