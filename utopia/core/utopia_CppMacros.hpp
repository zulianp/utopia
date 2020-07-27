#ifndef UTOPIA_CPP_MACROS_HPP
#define UTOPIA_CPP_MACROS_HPP
////////////////////////////////////////

#include <cassert>
#include "utopia_Base.hpp"

// Cuda should have its own assert(...) function
#define UTOPIA_DEVICE_ASSERT(expr) assert((expr))

#ifdef WITH_CPP17
/////////////////////////////
//////// C++17 //////////////
#define UTOPIA_IF_CONSTEXPR if constexpr
#define UTOPIA_CONSTEXPR constexpr
#define UTOPIA_DEVICE_ASSERT_CXX14(expr) UTOPIA_DEVICE_ASSERT((expr))

/////////////////////////////
#else
#ifdef WITH_CPP14
/////////////////////////////
//////// C++14 //////////////
#define UTOPIA_IF_CONSTEXPR if
#define UTOPIA_CONSTEXPR constexpr
#define UTOPIA_DEVICE_ASSERT_CXX14(expr) UTOPIA_DEVICE_ASSERT((expr))
/////////////////////////////

#else
/////////////////////////////
//////// C++11 //////////////
#define UTOPIA_IF_CONSTEXPR if
#define UTOPIA_CONSTEXPR
#define UTOPIA_DEVICE_ASSERT_CXX14(...)

/////////////////////////////
#endif  // WITH_CPP14
#endif  // WITH_CPP17

////////////////////////////////////////
#endif  // UTOPIA_CPP_MACROS_HPP
