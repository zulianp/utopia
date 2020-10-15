#include "utopia_AppRunner.hpp"
#include "utopia_AppRegistry.hpp"
#include "utopia_AppWithInputRegistry.hpp"
#include "utopia_Base.hpp"
#include "utopia_IOStream.hpp"
#include "utopia_Instance.hpp"
#include "utopia_NaryActionRegistry_impl.hpp"
#include "utopia_ui.hpp"

#include <iostream>

namespace utopia {

    template class NaryActionRegistry<Input &>;

    AppRunner::AppRunner() = default;

    AppRunner::~AppRunner() = default;

    void AppRunner::verbose(const bool val) { AppRegistry::instance().verbose(val); }

    int AppRunner::run(const std::string &name, Input &in) const {
        return AppWithInputRegistry::instance().apply(name, in);
    }

    int AppRunner::run(int argc, char **argv) {
        std::vector<std::string> apps;
        this->verbose(Utopia::instance().verbose());

        // utopia::out() <<argc << " " << argv[0] << " " << argv[1] << std::endl;

        int err = 0;

        if (argc == 2 && (argv[1] == std::string("-list") || argv[1] == std::string("-help"))) {
            this->describe();
            return 0;
        }

        if (argc == 4 && (argv[1] == std::string("-app"))) {
            auto in = open_istream(argv[3]);

            if (in) {
                if ((err = run(argv[2], *in)) == 0) {
                    std::cerr << "[Warning] syntax deprected use: " << argv[0] << " " << argv[1] << " " << argv[2]
                              << " @file " << argv[3] << std::endl;
                    return 0;
                }
            }
        }

        InputParameters params;
        params.init(argc, argv);

        std::string app_name = "";
        params.get("app", app_name);

        // try to run with parameters
        if ((err = run(app_name, params)) != 0) {
            // try to run without. parameters
            if ((err = run({app_name})) != 0) {
                std::cerr << "[Error] no app with name " << app_name << std::endl;
                return err;
            }
        }

        return err;
    }

    int AppRunner::run(const std::vector<std::string> &apps) const {
        auto &tr = AppRegistry::instance();
        int error_code = 0;
        for (const auto &t : apps) {
            int temp = 0;
            if ((temp = tr.run(t)) != 0) {
                error_code = temp;
            }
        }

        auto &awi = AppWithInputRegistry::instance();
        InputParameters empty;
        for (const auto &t : apps) {
            int temp = 0;
            if ((temp = awi.apply(t, empty)) != 0) {
                error_code = temp;
            }
        }

        return error_code;
    }

    void AppRunner::describe(std::ostream &os) const {
        utopia::out() << "Apps without input: " << std::endl;
        AppRegistry::instance().describe(os);

        utopia::out() << "Apps with input: " << std::endl;
        AppWithInputRegistry::instance().describe(os);
    }

}  // namespace utopia
