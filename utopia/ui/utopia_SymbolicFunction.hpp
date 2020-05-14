#ifndef UTOPIA_SYMBOLIC_FUNCTION_HPP
#define UTOPIA_SYMBOLIC_FUNCTION_HPP

#include "utopia_Base.hpp"

#ifdef WITH_TINY_EXPR

#include <memory>
#include "utopia_Expression.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    /**
     * @brief Symbolic evaluation of strings using the tinyexpr library
     * supported variable names are x, y, z, t and are ordered by dimension.
     */
    class SymbolicFunction final : public Expression<SymbolicFunction> {
    public:
        static const int Order = 0;
        using Scalar = double;

        class Impl;

        double eval(const double x);
        double eval(const double x, const double y);
        double eval(const double x, const double y, const double z);
        double eval(const double x, const double y, const double z, const double t);
        double eval(const std::vector<double> &x);

        bool valid() const;

        SymbolicFunction & operator=(const SymbolicFunction &other);
        SymbolicFunction(const std::string &expr);
        SymbolicFunction(const SymbolicFunction &other);
        ~SymbolicFunction();

    private:
        std::unique_ptr<Impl> impl_;
    };

    inline SymbolicFunction symbolic(const std::string &expr)
    {
        return SymbolicFunction(expr);
    }

    template<>
    class Traits<SymbolicFunction> : public Traits<double> {};

    template<class Left>
    class MostDescriptive<Left, SymbolicFunction > {
    public:
        using Type = Left;
    };

    template<class Right>
    class MostDescriptive<SymbolicFunction, Right > {
    public:
        using Type = Right;
    };

    template<>
    class MostDescriptive<SymbolicFunction, SymbolicFunction> {
    public:
        using Type = double;
    };

}

#endif //WITH_TINY_EXPR

#endif //UTOPIA_SYMBOLIC_FUNCTION_HPP
