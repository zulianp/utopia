#ifndef UTOPIA_CL_EVALUATOR_HPP
#define UTOPIA_CL_EVALUATOR_HPP 

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Operators.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Reduce.hpp"
#include "utopia_Backend.hpp"
#include "utopia_Assign.hpp"
#include "utopia_Structured.hpp"
#include "utopia_Factory.hpp"
#include "utopia_Ranged.hpp"
#include "utopia_Multiply.hpp"
#include "utopia_Transposed.hpp"
#include "utopia_Boolean.hpp"
#include "utopia_OuterProduct.hpp"
#include "utopia_Norm.hpp"

namespace utopia {

	namespace opencl {
		template<class Expr>
		class Program;

		template<class Expr>
		void eval_program(const Expr &expr)
		{
			Program<Expr>::instance().initialize(expr);
			Program<Expr>::instance().execute(expr);
		}
	}

	template<class Result>
	class Evaluator<Result, utopia::OPENCL_TAG> {
	public:
		typedef utopia::Traits<Result> Traits;

		typedef typename Traits::Vector Vector;
		typedef typename Traits::Matrix Matrix;
		typedef typename Traits::Scalar Scalar;

		template<class Derived>
		Scalar eval(const Expression<Derived> &expr)
		{   
			Number<Scalar> result(0);
			utopia::opencl::eval_program(
				Construct<Number<Scalar>, Derived>(result, expr.derived())
						);
			return result;
		}

		template<class Tensor>
		void eval(const Wrapper<Tensor, 1> &expr, Size &size)
		{
			size.set_dims(1);
		  	size.set(0, raw_type(expr).entries.size());
		}

		template<class Tensor>
		void eval(const Wrapper<Tensor, 2> &expr, Size &size)
		{
			size.set_dims(2);
		  	size.set(0, raw_type(expr).rows);
		  	size.set(1, raw_type(expr).cols);
		}

		template<class Left, class Right>
		void eval(const Assign< Wrapper<Left, 1>, Right> &expr)
		{
			Size s = size(expr.right());
			raw_type(expr.left()).resize(s.get(0));		
		  	utopia::opencl::eval_program(expr);  
		}

	
		template<class Left, class Right>
		void eval(const Construct< Wrapper<Left, 1>, Right> &expr)
		{
			Size s = size(expr.right());
			raw_type(expr.left()).resize(s.get(0));			
		    utopia::opencl::eval_program(expr);  
		}

		template<class Left, class Right>
		void eval(const Assign< Wrapper<Left, 2>, Right> &expr)
		{
			Size s = size(expr.right());
			raw_type(expr.left()).resize(s.get(0), s.get(1));		
		  	utopia::opencl::eval_program(expr);  
		}

		
		template<class Left, class Right>
		void eval(const Construct< Wrapper<Left, 2>, Right> &expr)
		{
			Size s = size(expr.right());
			raw_type(expr.left()).resize(s.get(0), s.get(1));			
		    utopia::opencl::eval_program(expr);  
		}
	};
}

#endif //UTOPIA_CL_EVALUATOR_HPP
