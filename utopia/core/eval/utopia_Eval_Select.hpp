#ifndef UTOPIA_EVAL_SELECT_HPP
#define UTOPIA_EVAL_SELECT_HPP

#include "utopia_Eval_Empty.hpp" 

namespace utopia {
	template<class Left, class Right, typename SizeType, class Traits, int Backend>
	class Eval< Assign<Left, Select<Right, SizeType, 1> >, Traits, Backend> {
	public:
	    inline static bool apply(const Assign<Left, Select<Right, SizeType, 1> > &expr)
	    {
	        UTOPIA_LOG_BEGIN(expr);

	        UTOPIA_BACKEND(Traits).select(
	                Eval<Left,  Traits>::apply(expr.left()),
	                Eval<Right, Traits>::apply(expr.right().expr()),
	                expr.right().index()
	        );

	        UTOPIA_LOG_END(expr);
	        return true;
	    }
	};


	template<class Left, class Right, typename SizeType, class Traits, int Backend>
	class Eval< Construct<Left, Select<Right, SizeType, 1> >, Traits, Backend> {
	public:
	    inline static bool apply(const Construct<Left, Select<Right, SizeType, 1> > &expr)
	    {
	        UTOPIA_LOG_BEGIN(expr);

	        UTOPIA_BACKEND(Traits).select(
	                Eval<Left,  Traits>::apply(expr.left()),
	                Eval<Right, Traits>::apply(expr.right().expr()),
	                expr.right().index()
	        );

	        UTOPIA_LOG_END(expr);
	        return true;
	    }
	};
}

#endif //UTOPIA_EVAL_SELECT_HPP
