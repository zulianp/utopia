#ifndef UTOPIA_TESTING_HPP
#define UTOPIA_TESTING_HPP

#include "utopia_TestMacros.hpp"
#include "utopia_TestRegistry.hpp"
#include "utopia_TestRunner.hpp"

#include "utopia_Chrono.hpp"
#include "utopia_Instance.hpp"
#include "utopia_MPI.hpp"

namespace utopia {
    int Test(const int argc, char *argv[]);
}  // namespace utopia

#define UTOPIA_TEST(macro_argc__, macro_argv__) utopia::Test(argc, argv)

#endif  // UTOPIA_TESTING_HPP
