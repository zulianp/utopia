#ifndef UTOPIA_BLAS_EVAL_MISC_HPP
#define UTOPIA_BLAS_EVAL_MISC_HPP

#include "utopia_Eval_Empty.hpp"

namespace utopia {

	//specialize this for blas
	template<class T, class Traits>
	class Eval< Binary<Transposed<T>, T, Plus>, Traits, utopia::BLAS > {
	public:
		using SizeType = typename Traits::SizeType;
	    using Expr = utopia::Binary<Transposed<T>, T, Plus>;
	    using Result = EXPR_TYPE(Traits, Expr);

	    inline static void apply(const Expr &expr, Result &result)
	    {
	    	UTOPIA_TRACE_BEGIN(expr);

	        auto &&left  = Eval<T, Traits>::apply(expr.left().expr());
	        auto &&right = Eval<T, Traits>::apply(expr.right());

	        assert(left.rows() == right.cols());
	        assert(left.cols() == right.rows());

	        if(result.is_alias(left)) {
	        	
	        	if(result.is_alias(right)) {
	        		result.add_transpose(result);
	        	} else {
	            	result.transpose_add(right);
	            }

	        } else {
	            apply_aux(left, right, result);
	        } 

	        UTOPIA_TRACE_END(expr);
	    }

	private:
	    template<class T1, class T2>
	    inline static void apply_aux(T1 &&left, T2 &&right, Result &result)
	    {
	        const auto rows = right.rows();
	        const auto cols = right.cols();

	        if(result.rows() != rows || result.cols() != cols) {
	            result.resize(rows, cols);
	        }

	        for(SizeType r = 0; r < rows; ++r) {
	            for(SizeType c = 0; c < cols; ++c) {
	                result.set(r, c, left.get(c, r) + right.get(r, c));
	            }
	        }
	    }

	};

	template<class T, class Traits>
	class Eval< Binary<T, Transposed<T>, Plus>, Traits, utopia::BLAS > {
	public:
	    using Expr = utopia::Binary<T, Transposed<T>, Plus>;
	    using Result = EXPR_TYPE(Traits, Expr);

	    inline static void apply(const Expr &expr, Result &result)
	    {
	        UTOPIA_TRACE_BEGIN(expr);

	        auto &&left  = Eval<T, Traits>::apply(expr.left());
	        auto &&right = Eval<T, Traits>::apply(expr.right().expr());

	        if(result.is_alias(left)) {
	            result.add_transpose(right);
	        } else if(result.is_alias(right)) {
	            result.transpose_add(left);
	        } else {
	            result.construct( left );
	            result.add_transpose(left);
	        }

	        UTOPIA_TRACE_END(expr);
	    }

	};

	template<class T, class Traits>
	class Eval< InPlace<T, Transposed<T>, Plus>, Traits, BLAS > {
	public:

		inline static void apply(const InPlace<T, Transposed<T>, Plus> &expr)
		{
			UTOPIA_TRACE_BEGIN(expr);

			auto &&left  = Eval<T, Traits>::apply(expr.left());
			auto &&right = Eval<T, Traits>::apply(expr.right().expr());

			assert(left.rows() == right.cols());
			assert(left.cols() == right.rows());

			const auto rows = left.rows();
			const auto cols = left.cols();

			left.add_transpose(right);

			UTOPIA_TRACE_END(expr);
		}
	};

}

#endif //UTOPIA_BLAS_EVAL_MISC_HPP
