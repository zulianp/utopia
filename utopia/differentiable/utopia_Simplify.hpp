#ifndef UTOPIA_SIMPLIFY_HPP
#define UTOPIA_SIMPLIFY_HPP

namespace utopia {

    template <class Expr>
    class Simplify {
    public:
        using Type = Expr;

        inline static const Expr &make(const Expr &expr) { return expr; }
    };
}  // namespace utopia

#endif  // UTOPIA_SIMPLIFY_HPP
