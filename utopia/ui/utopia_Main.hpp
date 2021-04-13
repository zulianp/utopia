#ifndef UTOPIA_MAIN_HPP
#define UTOPIA_MAIN_HPP

#include "utopia.hpp"
#include "utopia_AppRunner.hpp"
#include "utopia_Options.hpp"

namespace utopia {

    int Main(const int argc, char *argv[]);

}  // namespace utopia

#define UTOPIA_MAIN(macro_argc__, macro_argv__) utopia::Main(argc, argv)

#endif  // UTOPIA_MAIN_HPP
