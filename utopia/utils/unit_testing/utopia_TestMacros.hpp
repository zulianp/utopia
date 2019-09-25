#ifndef UTOPIA_TEST_MACROS_HPP
#define UTOPIA_TEST_MACROS_HPP

#include "utopia_AutoRegisterTestUnit.hpp"
#include "utopia_TestRegistry.hpp"

namespace utopia {

    class UnitTestBase {
    public:
        virtual ~UnitTestBase() {}
        virtual void set_up() {}
        virtual void tear_down() {}
    };

    #define UTOPIA_DEFINE_VAR(macro_in) dummy_test_variable_ ## macro_in ## __LINE__

    // #define UTOPIA_REGISTER_TEST_CLASS(TestClassName_) static utopia::AutoRegisterTestUnit<TestClassName_> dummy_test_variable_ ## TestClassName_ ## __LINE__(#TestClassName_)
    // #define UTOPIA_REGISTER_TEST_FUNCTION(TestCFunctionName_) static utopia::AutoRegisterTestUnitFunction dummy_test_variable_ ## TestCFunctionName_ ## __LINE__(#TestCFunctionName_, TestCFunctionName_)

    #define UTOPIA_REGISTER_TEST_FUNCTION(TestCFunctionName_) static char UTOPIA_DEFINE_VAR(TestCFunctionName_) = utopia::TestRegistry::add_test_unit(#TestCFunctionName_, TestCFunctionName_)
}


#endif //UTOPIA_TEST_MACROS_HPP
