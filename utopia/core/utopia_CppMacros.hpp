#ifndef UTOPIA_CPP_MACROS_HPP
#define UTOPIA_CPP_MACROS_HPP


#ifdef WITH_CPP17
/////////////////////////////
//////// C++17 //////////////

#define UTOPIA_IF_CONSTEXPR if constexpr

/////////////////////////////

#else
#ifdef WITH_CPP14
/////////////////////////////
//////// C++14 //////////////
#define UTOPIA_IF_CONSTEXPR if

/////////////////////////////

#else
/////////////////////////////
//////// C++11 //////////////
#define UTOPIA_IF_CONSTEXPR if

/////////////////////////////
#endif //WITH_CPP14
#endif //WITH_CPP17


#endif