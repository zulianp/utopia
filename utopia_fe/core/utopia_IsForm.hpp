#ifndef UTOPIA_IS_FORM_HPP
#define UTOPIA_IS_FORM_HPP

#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_FEIsSubTree.hpp"

namespace utopia {
    ///gives false positives
    template<class Expr>
    struct IsForm
    {
        typedef utopia::Reduce< Binary<utopia::Any, utopia::Any, EMultiplies>, Plus> Inner;
        typedef Binary<utopia::Any, utopia::Any, Multiplies> Inner2;

        static const int has_inner_product = IsSubTree<Inner, Expr>::value;
        static const int has_scalar_inner_product = IsSubTree<Inner2, Expr>::value;

        static const int has_trial = IsSubTree<TrialFunction<utopia::Any>, Expr>::value;
        static const int has_test  = IsSubTree<TestFunction<utopia::Any>,  Expr>::value;
        static const int has_fun   = has_trial || has_test;
        static const int order     = has_inner_product * (has_trial + has_test);
        static const int value     = has_inner_product && (has_fun);
    };


    template<class Left, class Right>
    struct IsBilinearFormSplit
    {
        //check if left has a fe function
        static const int left_has_trial = IsSubTree<TrialFunction<utopia::Any>, Left>::value;
        static const int left_has_test  = IsSubTree<TestFunction<utopia::Any>,  Left>::value;

        //check if right has a fe function
        static const int right_has_trial = IsSubTree<TrialFunction<utopia::Any>, Right>::value;
        static const int right_has_test  = IsSubTree<TestFunction<utopia::Any>,  Right>::value;
        static const int has_test = left_has_test | right_has_test;

        static const bool is_linear = (right_has_test || left_has_test) && (!left_has_trial && !right_has_trial);
        static const bool value = (left_has_trial && right_has_test) || (right_has_trial && left_has_test);
        static const int  order  = is_linear? 1 : (value? 2 : 0);


        static const int has_fun = left_has_trial || right_has_trial || left_has_test || right_has_test;
    };

    template<class Expr>
    struct IsBilinearForm
    {
        static const int value = IsForm<Expr>::order == 2;
    };

    template<class Expr>
    struct IsLinearForm
    {
        static const int value = IsForm<Expr>::order == 2;
    };


    template<class Expr>
    struct IsForm<const Expr> : IsForm<Expr> {};

    template<class Expr>
    struct IsForm<const Expr &> : IsForm<Expr> {};

    template<class Expr>
    struct IsForm<Expr &> : IsForm<Expr> {};

    template<class Expr>
    struct IsForm<Expr &&> : IsForm<Expr> {};
}

#endif //UTOPIA_IS_FORM_HPP
