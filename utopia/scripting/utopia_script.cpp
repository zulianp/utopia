#include "utopia_script.hpp"
#include "../utopia.hpp"
#include "utopia_Instance.hpp"
#include "utopia_Layout.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_Version.hpp"

#include <iostream>

namespace scripting {

    /////////////////////////////////Global////////////////////////////////////////////////

    void init(int argc, char *argv[]) {
        using namespace utopia;
        Utopia::Init(argc, argv);
        scripting::print_info();

        scripting::Factory::print_info();
    }

    void init() {
        // FIXME?
        int argc = 1;
        std::string argv = "utopia_script";
        char *argv_ptr = &argv[0];
        init(argc, &argv_ptr);
    }

    void print_info() { utopia::out() << "Utopia\nversion: " << UTOPIA_VERSION << std::endl; }

    void finalize() { utopia::Utopia::Finalize(); }

}  // namespace scripting