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
            bool ok = UTOPIA_BACKEND(Traits).transpose(
                    Eval<Tensor, Traits>::apply(t.expr()),
                    result);

            assert(ok);
            return result;
        }
    };
}

#endif //UTOPIA_UTOPIA_EVAL_TRANSPOSE_HPP
