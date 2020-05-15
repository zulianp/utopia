#ifndef UTOPIA_MACROS_HPP
#define UTOPIA_MACROS_HPP

// for avoiding rewriting the same code over and over again
#define UTOPIA_EVAL_APPLY_TO_TEMPORARY(Macro_Expr_, Macro_Result_) \
    inline static Macro_Result_ apply(const Macro_Expr_ &expr) {   \
        Macro_Result_ result;                                      \
        apply(expr, result);                                       \
        return result;                                             \
    }

#endif  // UTOPIA_MACROS_HPP
