#ifndef UTOPIA_APP_MACROS_HPP
#define UTOPIA_APP_MACROS_HPP

#include "utopia_AppRegistry.hpp"

#define UTOPIA_RUN_APP(app_name) \
    {                               \
        utopia::Chrono private_c; private_c.start(); \
        if(utopia::mpi_world_rank() == 0 && utopia::Utopia::instance().verbose()) { std::cout << "> " << std::left << std::setw(40) << (#app_name) << std::flush; } \
        app_name();                                \
        private_c.stop();                             \
        if(utopia::mpi_world_rank() == 0 && utopia::Utopia::instance().verbose()) { std::cout << "(" << private_c.get_seconds() << "s)" << std::endl; } \
    }

namespace utopia {

    #define UTOPIA_DEFINE_APP_VAR(macro_in) dummy_app_variable_ ## macro_in ## __LINE__
    #define UTOPIA_REGISTER_APP(AppCFunctionName_) static char UTOPIA_DEFINE_APP_VAR(AppCFunctionName_) = utopia::AppRunner::register_app(#AppCFunctionName_, AppCFunctionName_)
}

#endif //UTOPIA_APP_MACROS_HPP
