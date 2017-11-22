#ifndef UTOPIA_LIBMESH_LINEAR_FORM_EVAL_HPP
#define UTOPIA_LIBMESH_LINEAR_FORM_EVAL_HPP

#include "utopia_Base.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"

#include "utopia_FormEval.hpp"
#include "utopia_FEEval.hpp"

#include "utopia_libmesh_AssemblyContext.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"


namespace utopia {
	template<class Expr>
	class OffsetFinder {
	public:
		template<class Any>
		inline constexpr static int visit(const Any &) { return TRAVERSE_CONTINUE; }
		
		template<class T>
		inline int visit(const TestFunction<T> &expr)
		{	
			//TODO
			return TRAVERSE_CONTINUE;
		}

		template<class T>
		inline int visit(const TestFunction<ProductFunctionSpace<T>> &expr)
		{
			//TODO
			return TRAVERSE_CONTINUE;
		}

		template<class T>
		inline int visit(const TrialFunction<T> &expr)
		{
			//TODO
			return TRAVERSE_CONTINUE;
		}

		template<class T>
		inline int visit(const TrialFunction<ProductFunctionSpace<T>> &expr)
		{
			//TODO
			return TRAVERSE_CONTINUE;
		}

		template<class ExprTree>
		inline void apply(const ExprTree &expr)
		{
			traverse(expr, *this);
		}

		OffsetFinder(AssemblyContext<LIBMESH_TAG> &ctx)
		: ctx(ctx), row_offset(0), col_offset(0)
		{}

		AssemblyContext<LIBMESH_TAG> &ctx;

		int row_offset;
		int col_offset;
	};

	template<class Form>
	class FormEval<Form, LIBMESH_TAG> {
	public:
		typedef utopia::Traits<LibMeshFunctionSpace> Traits;
		FormEval() { }

		template<class Expr, class Tensor>
		static void apply(
					const Integral<Expr> &expr, 
					Tensor &t, 
					AssemblyContext<LIBMESH_TAG> &ctx)
		{
			if(expr.has_block_id() && ctx.block_id() != expr.block_id()) {
				return;
			}

			auto &&r = FEEval<Integral<Expr>, Traits, LIBMESH_TAG>::apply(expr, ctx);
			t = r;
		}


		template<class Left, class Right, class Matrix, class Vector>
		static void apply(
					const Equality<Left, Right> &expr, 
					Wrapper<Matrix, 2> &mat,
					Wrapper<Vector, 1> &vec, 
					AssemblyContext<LIBMESH_TAG> &ctx)
		{
			apply(expr.left(),  mat, ctx);
			apply(expr.right(), vec, ctx);
		}


		template<class... Eq, class Matrix, class Vector>
		static void apply(
			const Equations<Eq...> &eqs,
			Wrapper<Matrix, 2> &mat,
			Wrapper<Vector, 1> &vec, 
			AssemblyContext<LIBMESH_TAG> &ctx)
		{
			FEBackend<LIBMESH_TAG>::init_context(eqs, ctx);
		}

		template<class Left, class Right, class Tensor>
		static void apply(
			const Binary<Left, Right, Plus> &expr, 
			Tensor &result, 
			AssemblyContext<LIBMESH_TAG> &ctx)
		{	
			apply(expr.left(),  result, ctx);

			Tensor right = zeros(size(result));
			apply(expr.right(), right, ctx);
			result += right;
		}

		template<class Left, class Right, class Tensor>
		static void apply(
			const Binary<Number<Left>, Right, Multiplies> &expr, 
			Tensor &result, 
			AssemblyContext<LIBMESH_TAG> &ctx)
		{	
			apply(expr.right(), result, ctx);
			result *= expr.left();
		}

		template<class Left, class Right, class Tensor>
		static void apply(
			const Binary<Left, Right, Minus> &expr, 
			Tensor &result, 
			AssemblyContext<LIBMESH_TAG> &ctx)
		{	
			apply(expr.left(), result, ctx);

			Tensor right = zeros(size(result));
			apply(expr.right(), right, ctx);	
			result -= right;
		}

		template<class Left, class Right, class Op, class Tensor>
		static void apply(
			const Binary<Left, Right, Op> &expr, 
			Tensor &result, 
			AssemblyContext<LIBMESH_TAG> &ctx)
		{	
			apply(expr.left(), result, ctx);

			Tensor right = zeros(size(result));
			apply(expr.right(), result, ctx);	

			result = expr.operation().apply(result, right);
		}

		template<class Expr, class Tensor>
		static void apply(
			const Unary<Expr, Minus> &expr, 
			Tensor &result, 
			AssemblyContext<LIBMESH_TAG> &ctx)
		{	
			apply(expr.expr(), result, ctx);
			result = -result;
		}

		template<class Expr, class Tensor>
		static void apply(
			const Negate<Expr> &expr, 
			Tensor &result, 
			AssemblyContext<LIBMESH_TAG> &ctx)
		{	
			apply(expr.expr(), result, ctx);
			result = -result;
		}

		template<class Expr, class Tensor>
		static void apply(
			const Unary<Expr, Abs> &expr, 
			Tensor &result, 
			AssemblyContext<LIBMESH_TAG> &ctx)
		{	
			apply(expr.expr(), result, ctx);
			result = abs(result);
		}

		template<class Expr, class Tensor>
		static void apply(
			const Unary<Expr, Sqrt> &expr, 
			Tensor &result, 
			AssemblyContext<LIBMESH_TAG> &ctx)
		{	
			apply(expr.expr(), result, ctx);
			result = sqrt(result);
		}
	};
}

#endif //UTOPIA_LIBMESH_LINEAR_FORM_EVAL_HPP
