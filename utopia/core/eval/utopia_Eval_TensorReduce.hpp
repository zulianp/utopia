#ifndef UTOPIA_EVAL_TENSOR_REDUCE_HPP
#define UTOPIA_EVAL_TENSOR_REDUCE_HPP

#include "utopia_Eval_Empty.hpp"
// #include "utopia_TensorReduce.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

	template<class Expr, class Operation, class Traits, int Backend>
	class Eval< TensorReduce<Expr, Operation>, Traits, Backend> {
	public:
	    typedef typename TypeAndFill<Traits,TensorReduce<Expr, Operation> >::Type Result;

	    inline static Result apply(const TensorReduce<Expr, Operation> &expr) {
	        Result result;

            UTOPIA_LOG_BEGIN(expr);

	        const bool ok = UTOPIA_BACKEND(Traits).apply_tensor_reduce(
	                Eval<Expr,  Traits>::apply(expr.expr()),
	                expr.operation(),
	                expr.dim(),
	                result
	        );

	        ASSERT(ok);

            UTOPIA_LOG_END(expr);
	        return result;
	    }
	};
}


#endif //UTOPIA_EVAL_TENSOR_REDUCE_HPP
