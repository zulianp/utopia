//
// Created by Patrick Zulian on 15/05/15.
//

#ifndef SIMMOD_utopia_NEGATE_HPP
#define SIMMOD_utopia_NEGATE_HPP

#include "utopia_Operators.hpp"
#include "utopia_Unary.hpp"

namespace utopia {
    template<class Expr>
    class Negate : public Unary<Expr, Minus> {
    public:
        Negate(const Expr &expr) : Unary<Expr, Minus>(expr) { }

        std::string getClass() const override 
        {
        	return "Negate<" + this->expr().getClass() + ">";
        }

        virtual ~Negate() { }
    };
}
#endif //SIMMOD_utopia_NEGATE_HPP
