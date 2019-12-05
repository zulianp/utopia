#include "utopia_AppRunner.hpp"
#include "utopia_AppRegistry.hpp"
#include "utopia_Base.hpp"

#include <iostream>

namespace utopia {

    AppRunner::AppRunner()
    {}

    AppRunner::~AppRunner()
    {}

    void AppRunner::verbose(const bool val)
    {
        AppRegistry::instance().verbose(val);
    }

    int AppRunner::run(int argc, char **argv) const {
        UTOPIA_UNUSED(argc);
        UTOPIA_UNUSED(argv);
        return AppRegistry::instance().run_all();
    }

    int AppRunner::run(const std::vector<std::string> &apps) const
    {
        auto &tr = AppRegistry::instance();
        int error_code = 0;
        for(const auto &t : apps) {
            int temp = 0;
            if((temp = tr.run(t)) != 0) {
                error_code = temp;
            }
        }

        return error_code;
    }

    void AppRunner::describe(std::ostream &os) const
    {
        AppRegistry::instance().describe(os);
    }

}

