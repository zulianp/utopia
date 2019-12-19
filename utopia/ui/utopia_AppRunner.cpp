#include "utopia_AppRunner.hpp"
#include "utopia_AppRegistry.hpp"
#include "utopia_AppWithInputRegistry.hpp"
#include "utopia_NaryActionRegistry_impl.hpp"
#include "utopia_Base.hpp"
#include "utopia_Instance.hpp"
#include "utopia_ui.hpp"

#include <iostream>

namespace utopia {

    template class NaryActionRegistry<Input &>;

    AppRunner::AppRunner()
    {}

    AppRunner::~AppRunner()
    {}

    void AppRunner::verbose(const bool val)
    {
        AppRegistry::instance().verbose(val);
    }

    int AppRunner::run(const std::string &name, Input &in) const
    {
        return AppWithInputRegistry::instance().apply(name, in);
    }

    int AppRunner::run(int argc, char **argv) {

        std::vector<std::string> apps;
        this->verbose(Utopia::instance().verbose());

        bool app_ran = false;
        for(int i = 1; i < argc; i++) {
          if(argv[i] == std::string("-app")) {
                if(++i >= argc)
                    break;

                int err = 1;
                if(i + 1 < argc) {
                    auto istr = open_istream(argv[i+1]);
                    if(istr) {
                       err = run(argv[i], *istr);
                       if(err == 0) { app_ran = true; }
                    } else {
                        err = 1;
                    }
                }

                //Try without input
                if(err != 0) {
                    apps.push_back(argv[i]);
                }

            } else if(argv[i] == std::string("-list")) {
                this->describe();
                return 0;
            }
        }


        if(!app_ran) {
            if(apps.empty()) {
                return AppRegistry::instance().run_all();
            } else {
                return this->run(apps);
            }
        } else {
            return 0;
        }
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

        auto &awi = AppWithInputRegistry::instance();
        InputParameters empty;
        for(const auto &t : apps) {
            int temp = 0;
            if((temp = awi.apply(t, empty)) != 0) {
                error_code = temp;
            }
        }

        return error_code;
    }

    void AppRunner::describe(std::ostream &os) const
    {
        std::cout << "Apps without input: " << std::endl;
        AppRegistry::instance().describe(os);

        std::cout << "Apps with input: " << std::endl;
        AppWithInputRegistry::instance().describe(os);
    }

}

