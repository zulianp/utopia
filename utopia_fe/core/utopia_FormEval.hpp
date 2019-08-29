#ifndef UTOPIA_FORM_EVAL_HPP
#define UTOPIA_FORM_EVAL_HPP

#include "utopia_Traits.hpp"
#include <iostream>

namespace utopia {
    template<class Form, int BAKEND_FLAG>
    class FormEval {
    public:
        FormEval()
        {
            static_assert(BAKEND_FLAG < utopia::HOMEMADE, "FormEval: unimplemented evaluator for backend");
        }

        template<class Expr, class Tensor, class Context>
        static void apply(const Expr &expr, Tensor &, const Context &) {
            std::cerr << "[Error] not implemented for expression " << std::endl;
            std::cerr << tree_format(expr.getClass()) << std::endl;
        }
    };

}

#endif //UTOPIA_FORM_EVAL_HPP
