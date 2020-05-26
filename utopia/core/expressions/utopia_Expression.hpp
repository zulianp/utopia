//
// Created by Patrick Zulian on 15/05/15.
//

#ifndef utopia_utopia_EXPRESSION_HPP
#define utopia_utopia_EXPRESSION_HPP

#include <string>
#include "utopia_StoreAs.hpp"

namespace utopia {

    /*!
     * @brief Every result of an algebraic
     * operation is an expression.
     * @tparam Derived is the most specific subclass of Expression.
     * It is using the Curiously Recurring Template (CRT) pattern.
     */
    template <class Derived>
    class Expression {
    public:
        enum { StoreAs = UTOPIA_DEFAULT_EXPRESSION_STORAGE };

        Derived &derived() { return static_cast<Derived &>(*this); }

        // virtual ~Expression() { }

        ///@return itself as its Derived type
        inline constexpr const Derived &derived() const { return static_cast<const Derived &>(*this); }

        ///@return itself as its Derived type
        inline operator Derived &() { return derived(); }

        ///@return itself as its Derived type
        inline constexpr operator const Derived &() const { return derived(); }

        ///@return a string with the name of the class
        virtual std::string get_class() const { return "Expression"; }
    };
}  // namespace utopia

#endif  // utopia_utopia_EXPRESSION_HPP