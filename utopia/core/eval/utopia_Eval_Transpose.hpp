//
// Created by Patrick Zulian on 29/08/16.
//

#ifndef UTOPIA_UTOPIA_EVAL_TRANSPOSE_HPP
#define UTOPIA_UTOPIA_EVAL_TRANSPOSE_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

    template<class Tensor, class Traits, int Backend>
    class Eval<Transposed<Tensor>, Traits, Backend> {
    public:
        typedef typename TypeAndFill<Traits, Tensor>::Type Result;

        inline static Result apply(const Transposed<Tensor> &t)
        {
            Result result;

            UTOPIA_TRACE_BEGIN(t);

            UTOPIA_BACKEND(Traits).assign_transposed(
                result,
                Eval<Tensor, Traits>::apply(t.expr())
                );

            UTOPIA_TRACE_END(t);
            return result;
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_TRANSPOSE_HPP
