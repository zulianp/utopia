#ifndef UTOPIA_EVAL_SELECT_HPP
#define UTOPIA_EVAL_SELECT_HPP

#include "utopia_Eval_Empty.hpp" 

namespace utopia {
	template<class Left, class Right, typename SizeType, class Traits, int Backend>
	class Eval< Assign<Left, Select<Right, SizeType, 1> >, Traits, Backend> {
	public:
	    inline static bool apply(const Assign<Left, Select<Right, SizeType, 1> > &expr)
	    {
	        UTOPIA_TRACE_BEGIN(expr);

	        UTOPIA_BACKEND(Traits).select(
	                Eval<Left,  Traits>::apply(expr.left()),
	                Eval<Right, Traits>::apply(expr.right().expr()),
	                expr.right().index()
	        );

	        UTOPIA_TRACE_END(expr);
	        return true;
	    }
	};


	template<class Left, class Right, typename SizeType, class Traits, int Backend>
	class Eval< Construct<Left, Select<Right, SizeType, 1> >, Traits, Backend> {
	public:
	    inline static bool apply(const Construct<Left, Select<Right, SizeType, 1> > &expr)
	    {
	        UTOPIA_TRACE_BEGIN(expr);

	        UTOPIA_BACKEND(Traits).select(
	                Eval<Left,  Traits>::apply(expr.left()),
	                Eval<Right, Traits>::apply(expr.right().expr()),
	                expr.right().index()
	        );

	        UTOPIA_TRACE_END(expr);
	        return true;
	    }
	};

	template<class Expr, typename SizeType, class Traits, int Backend>
	class Eval< Select<Expr, SizeType, 1>, Traits, Backend> {
	public:
		typedef typename TypeAndFill<Traits, Expr>::Type Result;

	    inline static Result apply(const Select<Expr, SizeType, 1> &expr)
	    {
	       UTOPIA_TRACE_BEGIN(expr);
	       Result result;
	       
	       UTOPIA_BACKEND(Traits).select(
	               result,
	               Eval<Expr, Traits>::apply(expr.expr()),
	               expr.index()
	       );

	        UTOPIA_TRACE_END(expr);
	        return result;
	    }
	};


	template<class Left, class Right, typename SizeType, class Traits, int Backend>
	class Eval< Assign<Left, Select<Right, SizeType, 2> >, Traits, Backend> {
	public:
	    inline static bool apply(const Assign<Left, Select<Right, SizeType, 2> > &expr)
	    {
	        UTOPIA_TRACE_BEGIN(expr);

	        UTOPIA_BACKEND(Traits).select(
	                Eval<Left,  Traits>::apply(expr.left()),
	                Eval<Right, Traits>::apply(expr.right().expr()),
	                expr.right().row_index(),
	                expr.right().col_index()
	        );

	        UTOPIA_TRACE_END(expr);
	        return true;
	    }
	};


	template<class Left, class Right, typename SizeType, class Traits, int Backend>
	class Eval< Construct<Left, Select<Right, SizeType, 2> >, Traits, Backend> {
	public:
	    inline static bool apply(const Construct<Left, Select<Right, SizeType, 2> > &expr)
	    {
	        UTOPIA_TRACE_BEGIN(expr);

	       UTOPIA_BACKEND(Traits).select(
	               Eval<Left,  Traits>::apply(expr.left()),
	               Eval<Right, Traits>::apply(expr.right().expr()),
	               expr.right().row_index(),
	               expr.right().col_index()
	       );

	        UTOPIA_TRACE_END(expr);
	        return true;
	    }
	};


	template<class Expr, typename SizeType, class Traits, int Backend>
	class Eval< Select<Expr, SizeType, 2>, Traits, Backend> {
	public:
		typedef typename TypeAndFill<Traits, Expr>::Type Result;

	    inline static Result apply(const Select<Expr, SizeType, 2> &expr)
	    {
	       UTOPIA_TRACE_BEGIN(expr);
	       Result result;

	       UTOPIA_BACKEND(Traits).select(
	               result,
	               Eval<Expr, Traits>::apply(expr.expr()),
	               expr.row_index(),
	               expr.col_index()
	       );

	        UTOPIA_TRACE_END(expr);
	        return result;
	    }
	};
}

#endif //UTOPIA_EVAL_SELECT_HPP
