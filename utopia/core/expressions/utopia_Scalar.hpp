//
// Created by Patrick Zulian on 15/05/15.
//

#ifndef utopia_utopia_SCALAR_HPP
#define utopia_utopia_SCALAR_HPP

#include "utopia_Expression.hpp"

namespace utopia {
    //    class Scalar : public Expression<Scalar> {
    //    public:
    //        std::string get_class() const { return "Scalar"; }
    //    };

    template <typename T>
    class Value : public Expression<Value<T> > {
    public:
        inline Value(const T &value) : _value(value) {}
        inline operator T() const { return _value; }

    private:
        T _value;
    };
}  // namespace utopia
#endif  // utopia_utopia_SCALAR_HPP
