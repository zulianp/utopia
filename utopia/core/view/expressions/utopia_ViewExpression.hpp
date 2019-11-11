#ifndef UTOPIA_VIEW_EXPRESSION_HPP
#define UTOPIA_VIEW_EXPRESSION_HPP

#include "utopia_Base.hpp"
#include "utopia_StoreAs.hpp"
#include <string>

namespace utopia {

    /*!
     * @brief Every result of an algebraic
     * operation is an expression.
     * @tparam Derived is the most specific subclass of ViewExpression.
     * It is using the Curiously Recurring Template (CRT) pattern.
     */
    template<class Derived>
    class ViewExpression {
    public:
        enum {
            StoreAs = UTOPIA_DEFAULT_EXPRESSION_STORAGE
        };

        UTOPIA_INLINE_FUNCTION Derived &derived() { return static_cast<Derived &>(*this); }

        // virtual ~ViewExpression() { }

        ///@return itself as its Derived type
        UTOPIA_INLINE_FUNCTION constexpr const Derived &derived() const {
            return static_cast<const Derived &>(*this);
        }

        ///@return itself as its Derived type
        UTOPIA_INLINE_FUNCTION operator Derived &() {
            return derived();
        }

        ///@return itself as its Derived type
        UTOPIA_INLINE_FUNCTION constexpr  operator const Derived &() const {
            return derived();
        }

        ///@return a string with the name of the class
        virtual std::string get_class() const {
            return "ViewExpression";
        }
    };
}

#endif //UTOPIA_VIEW_EXPRESSION_HPP