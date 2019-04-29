#ifndef UTOPIA_SIMPLIFY_HPP
#define UTOPIA_SIMPLIFY_HPP

namespace utopia {

    template<class Expr>
    class Simplify {
    public:
        typedef Expr Type;

        inline static const Expr &make(const Expr &expr)
        {
            return expr;
        }
    };
}

#endif //UTOPIA_SIMPLIFY_HPP

