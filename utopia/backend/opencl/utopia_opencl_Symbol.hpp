#ifndef UTOPIA_SYMBOL_HPP
#define UTOPIA_SYMBOL_HPP

#include "utopia_Operators.hpp"

namespace utopia {
    namespace opencl {
        template<class Operation>
        class Symbol{};

        template<>
        class Symbol<utopia::Plus> {
        public:
            constexpr static const char *fun_str()
            {
                return "plus";
            }

            constexpr static const char *str()
            {
                return "+";
            }
        };

        template<>
        class Symbol<utopia::Minus> {
        public:
            constexpr static const char *fun_str()
            {
                return "minus";
            }


            constexpr static const char *str()
            {
                return "-";
            }
        };

        template<>
        class Symbol<utopia::Multiplies> {
        public:
            constexpr static const char *fun_str()
            {
                return "multiplies";
            }

            constexpr static const char *str()
            {
                return "*";
            }
        };


        template<>
        class Symbol<utopia::EMultiplies> {
        public:
            constexpr static const char *fun_str()
            {
                return "multiplies";
            }

            constexpr static const char *str()
            {
                return "*";
            }
        };



        template<>
        class Symbol<utopia::Divides> {
        public:
            constexpr static const char *fun_str()
            {
                return "divides";
            }

            constexpr static const char *str()
            {
                return "/";
            }
        };

        template<>
        class Symbol<utopia::Sqrt> {
        public:
            constexpr static const char *fun_str()
            {
                return "sqrt";
            }

            constexpr static const char *str()
            {
                return "sqrt";
            }
        };

        template<>
        class Symbol<utopia::Pow2> {
        public:
            constexpr static const char *fun_str()
            {
                return "pow2";
            }

            constexpr static const char *str()
            {
                return "pow2";
            }
        };

        template<>
        class Symbol<utopia::Abs> {
        public:
            constexpr static const char *fun_str()
            {
                return "fabs";
            }

            constexpr static const char *str()
            {
                return "fabs";
            }
        };
    }
}

#endif //UTOPIA_SYMBOL_HPP
