#ifndef UTOPIA_TYPE_TO_STRING_HPP
#define UTOPIA_TYPE_TO_STRING_HPP

#include "utopia_Path.hpp"

namespace utopia {
    template <typename T>
    class TypeToString {
    public:
        static constexpr const char *get() { return "object"; }
    };

#define UTOPIA_DEFINE_SUBS_TYPE_TO_STRING(T, Name)          \
    template <>                                             \
    class TypeToString<T> {                                 \
    public:                                                 \
        static constexpr const char *get() { return Name; } \
    };

#define UTOPIA_DEFINE_TYPE_TO_STRING(T) UTOPIA_DEFINE_SUBS_TYPE_TO_STRING(T, #T)

    UTOPIA_DEFINE_TYPE_TO_STRING(double);
    UTOPIA_DEFINE_TYPE_TO_STRING(float);
    UTOPIA_DEFINE_TYPE_TO_STRING(int);
    UTOPIA_DEFINE_TYPE_TO_STRING(bool);
    UTOPIA_DEFINE_TYPE_TO_STRING(long);
    UTOPIA_DEFINE_TYPE_TO_STRING(unsigned int);
    UTOPIA_DEFINE_TYPE_TO_STRING(unsigned long);
    UTOPIA_DEFINE_TYPE_TO_STRING(char);

    UTOPIA_DEFINE_SUBS_TYPE_TO_STRING(std::string, "string");
    UTOPIA_DEFINE_SUBS_TYPE_TO_STRING(utopia::Path, "path");

// Clean-up macros
#undef UTOPIA_DEFINE_SUBS_TYPE_TO_STRING
#undef UTOPIA_DEFINE_TYPE_TO_STRING

}  // namespace utopia

#endif