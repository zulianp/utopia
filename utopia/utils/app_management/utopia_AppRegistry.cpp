#include "utopia_AppRegistry.hpp"
#include "utopia_Chrono.hpp"
#include "utopia_MPI.hpp"

#include <iostream>

namespace utopia {

    AppRegistry::AppRegistry() { apps_.set_type("app"); }

    AppRegistry::~AppRegistry() = default;

    bool AppRegistry::verbose() const { return apps_.verbose(); }
    void AppRegistry::verbose(const bool val) { apps_.verbose(val); }

    char AppRegistry::add_app(const std::string &app_name, AppRegistry::RunApp run_app) {
        instance().apps_.add_action(app_name, run_app);
        return 0;
    }

    AppRegistry &AppRegistry::instance() {
        static AppRegistry instance_;
        return instance_;
    }

    void AppRegistry::describe(std::ostream &os) const {
        os << "Number of apps units: " << apps_.size() << std::endl;
        os << "select with: -app <sub-command>\n";
        os << "available sub-commands:\n";
        apps_.describe(os);
        os << std::flush;
    }

    int AppRegistry::run_all() { return apps_.apply_all(); }

    int AppRegistry::run(const std::string &app_name) {
        int ret = apps_.apply(app_name);
        // if(ret == -1) {
        //     std::cerr << "[Error] no app with name " << app_name << std::endl;
        // }

        return ret;
    }

}  // namespace utopia
