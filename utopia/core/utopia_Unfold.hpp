#ifndef UTOPIA_UNFOLD_HPP
#define UTOPIA_UNFOLD_HPP

namespace utopia {
    template <class Expr>
    class Unfold {
    public:
        using Type = Expr;
    };

    template <class Expr>
    class Fold {
        using Type = Expr;
    };
}  // namespace utopia

#endif  // UTOPIA_UNFOLD_HPP
