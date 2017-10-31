#ifndef UTOPIA_HOMEMADE_LINEAR_FORM_EVAL_HPP
#define UTOPIA_HOMEMADE_LINEAR_FORM_EVAL_HPP

#include "utopia_Base.hpp"
#include "utopia_homemade_FEForwardDeclarations.hpp"

#include "utopia_FormEval.hpp"
#include "utopia_FEEval.hpp"

#include "utopia_homemade_AssemblyContext.hpp"
#include "utopia_homemade_FunctionSpace.hpp"


namespace utopia {
	template<class Form>
	class FormEval<Form, HOMEMADE> {
	public:
		typedef utopia::Traits<HMFESpace> Traits;
		FormEval() { }

		template<typename T>
		inline static const T &get(const std::vector<std::vector<T> > &v, const std::size_t qp, const std::size_t i)
		{
			return v[qp][i];
		}

		template<typename T, int Order>
		inline static const T get(const ConstantCoefficient<T, Order> &c, const std::size_t qp, const std::size_t i)
		{
			return c[i];
		}

		template<typename T>
		inline static void add(ElementMatrix &mat, const int i, const int j, const T value)
		{
			mat.add(i, j, value);
		}

		template<typename T>
		inline static void add(ElementVector &vec, const int i, const int j, const T value)
		{
			assert(j == 0);
			vec.add(i, value);
		}

			/////////////////////////////////
			//FIXME extract to "super-class"
		template<class Left, class Right, class Tensor>
		static void apply(
			const Binary<Left, Right, Plus> &expr, 
			Tensor &result, 
			AssemblyContext<HOMEMADE> &ctx)
		{	
			apply(expr.left(),  result, ctx);

			Tensor right = zeros(size(result));
			apply(expr.right(), result, ctx);
			result += right;
		}

		template<class Left, class Right, class Tensor>
		static void apply(
			const Binary<Number<Left>, Right, Multiplies> &expr, 
			Tensor &result, 
			AssemblyContext<HOMEMADE> &ctx)
		{	
			apply(expr.right(), result, ctx);
			result *= expr.left();
		}

		template<class Left, class Right, class Tensor>
		static void apply(
			const Binary<Left, Right, Minus> &expr, 
			Tensor &result, 
			AssemblyContext<HOMEMADE> &ctx)
		{	
			apply(expr.left(), result, ctx);

			Tensor right = zeros(size(result));
			apply(expr.right(), result, ctx);	
			result -= right;
		}

		template<class Left, class Right, class Op, class Tensor>
		static void apply(
			const Binary<Left, Right, Op> &expr, 
			Tensor &result, 
			AssemblyContext<HOMEMADE> &ctx)
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
			AssemblyContext<HOMEMADE> &ctx)
		{	
			apply(expr.expr(), result, ctx);
			result = -result;
		}

		template<class Expr, class Tensor>
		static void apply(
			const Unary<Expr, Abs> &expr, 
			Tensor &result, 
			AssemblyContext<HOMEMADE> &ctx)
		{	
			apply(expr.expr(), result, ctx);
			result = abs(result);
		}

		template<class Expr, class Tensor>
		static void apply(
			const Unary<Expr, Sqrt> &expr, 
			Tensor &result, 
			AssemblyContext<HOMEMADE> &ctx)
		{	
			apply(expr.expr(), result, ctx);
			result = sqrt(result);
		}

			/////////////////////////////////


		template<class Expr, class Tensor>
		static void apply(
			const Integral<Expr> &expr, 
			Tensor &mat, 
			AssemblyContext<HOMEMADE> &ctx)
		{
			if(expr.has_block_id() && ctx.block_id() != expr.block_id()) {
				return;
			}

			apply(expr.expr(), mat, ctx);
		}



		template<class Left, class Right, class Tensor>
		static void apply(
			const Reduce<Binary<Left, Right, EMultiplies>, Plus> &expr, 
			Tensor &result, 
			AssemblyContext<HOMEMADE> &ctx)
		{	
			Write<Tensor> wt(result);

			auto && left  = FEEval<Left,  Traits, HOMEMADE>::apply(expr.expr().left(),  ctx);
			auto && right = FEEval<Right, Traits, HOMEMADE>::apply(expr.expr().right(), ctx);
			auto && dx    = ctx.dx();

			const bool left_is_test = is_test(expr.expr().left());
			assert( left_is_test != is_test(expr.expr().right()) );

			uint n_quad_points = dx.size();

			auto s = size(result);

			if(s.n_dims() == 1) {
				s.set_dims(2);
				s.set(1, 1);
			}

			if(left_is_test) {
				for (uint qp = 0; qp < n_quad_points; qp++) {
					for (uint i = 0; i < s.get(0); i++) {
						for (uint j = 0; j < s.get(1); j++) {
							add(result, i, j, inner( get(left, qp, i), get(right, qp, j) ) * dx[qp]);
						}
					}
				}

			} else {
				for (uint qp = 0; qp < n_quad_points; qp++) {
					for (uint i = 0; i < s.get(1); i++) {
						for (uint j = 0; j < s.get(0); j++) {
							add(result, j, i, inner( get(left, qp, i), get(right, qp, j) ) * dx[qp]);
						}
					}
				}
			}
		}

			///bilinear functional
		template<class Left, class Right>
		static Matrixd apply_bilinear(
			const Reduce<Binary<Left, Right, EMultiplies>, Plus> &expr, 
			AssemblyContext<HOMEMADE> &ctx)
		{
			Matrixd result;
			ctx.init_tensor(expr, result, true);
			apply(expr, result, ctx);
			return result;
		}


			///linear functional
		template<class Left, class Right>
		static Matrixd apply_linear(
			const Reduce<Binary<Left, Right, EMultiplies>, Plus> &expr, 
			AssemblyContext<HOMEMADE> &ctx)
		{
			Vectord result;
			ctx.init_tensor(expr, result, true);
			apply(expr, result, ctx);
			return result;
		}
	};
}

#endif //UTOPIA_HOMEMADE_LINEAR_FORM_EVAL_HPP
