// #ifndef UTOPIA_AUTO_REGISTER_TEST_UNIT_HPP
// #define UTOPIA_AUTO_REGISTER_TEST_UNIT_HPP

// #include "utopia_TestRegistry.hpp"
// #include <functional>

// namespace utopia {
    
//     template<class TestUnit>
//     class AutoRegisterTestUnit {
//     public:
//         AutoRegisterTestUnit(const char * name)
//         {
//             utopia::TestRegistry::instance().add_test_unit(name, &TestUnit::run);
//         }
//     };


//     class AutoRegisterTestUnitFunction {
//     public:
//         typedef void (*RunTest)();
//         AutoRegisterTestUnitFunction(const char * name, RunTest run)
//         {
//             utopia::TestRegistry::instance().add_test_unit(name, run);
//         }
//     };

// }


// #endif //UTOPIA_AUTO_REGISTER_TEST_UNIT_HPP
