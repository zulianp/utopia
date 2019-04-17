#ifndef UTOPIA_UNFOLD_HPP
#define UTOPIA_UNFOLD_HPP

namespace utopia {
    template<class Expr>
    class Unfold {
    public:
        typedef Expr Type;
    };


    template<class Expr>
    class Fold {
        typedef Expr Type;
    };
}

#endif //UTOPIA_UNFOLD_HPP
