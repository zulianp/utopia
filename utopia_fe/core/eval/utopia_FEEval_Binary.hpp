#ifndef UTOPIA_FE_EVAL_BINARY_HPP
#define UTOPIA_FE_EVAL_BINARY_HPP

#include "utopia_FEEval_Empty.hpp"
#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_IsForm.hpp"

namespace utopia {

	template<class Expr,
			 class Traits,
			 int IsQuadData,
			 int FormOrder = IsForm<Expr>::order,
			 int IsFE = IsForm<Expr>::has_fun || IsSubTree< Interpolate<utopia::Any, utopia::Any>, Expr>::value
			 >
	class BinaryDelegate {};


	template<class Expr, class Traits, int IsQuadData>
	class BinaryDelegate<Expr, Traits, IsQuadData, 0, 0> {
	public:

		template<class AssemblyContextT>
		static auto apply(const Expr &expr, AssemblyContextT &) -> const Expr &
		{
			return expr;
		}
	};

	template<class Expr, class Traits, int IsQuadData>
	class BinaryDelegate<Expr, Traits, IsQuadData, 0, 1> {
	public:

		template<class Left, class Right, class Op>
		static auto apply(const Binary<Left, Right, Op> &expr, AssemblyContext<Traits::Backend> &ctx) -> decltype(
			 FEBackend<Traits::Backend>::apply_binary(
							FEEval<Left,  Traits, Traits::Backend, IsQuadData>::apply(expr.left(),  ctx),
							FEEval<Right, Traits, Traits::Backend, IsQuadData>::apply(expr.right(), ctx),
							expr.operation(),
							ctx)
			 )
		{
			return FEBackend<Traits::Backend>::apply_binary(
				FEEval<Left, Traits,  Traits::Backend, IsQuadData>::apply(expr.left(),  ctx),
				FEEval<Right, Traits, Traits::Backend, IsQuadData>::apply(expr.right(), ctx),
				expr.operation(),
				ctx);
		}
	};

	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	template<class Left, class Right, class Op, class Traits, int IsQuadData>
	class BinaryDelegate<Binary<Left, Right, Op>, Traits, IsQuadData, 1, 1> {
	public:
		typedef utopia::Binary<Left, Right, Op> BinExpr;

		static const int Backend = Traits::Backend;
		typedef typename Traits::Vector Vector;

		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> Vector
		{
			return FEBackend<Backend>::apply_binary(
				FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx),
				FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(), ctx),
				expr.operation(),
				ctx);
		}
	};


	template<class Left, class Right, class Traits, int IsQuadData>
	class BinaryDelegate<Binary<Left, Right, Plus>, Traits, IsQuadData, 1, 1> {
	public:
		typedef utopia::Binary<Left, Right, Plus> BinExpr;

		static const int Backend = Traits::Backend;
		typedef typename Traits::Vector Vector;

		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> Vector
		{
			return FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx)  +
	   			   FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(), ctx);
		}
	};

	template<class Left, class Right, class Traits>
	class BinaryDelegate<Binary<Left, Right, Plus>, Traits, 1, 1, 1> {
	public:
		typedef utopia::Binary<Left, Right, Plus> BinExpr;

		static const int Backend = Traits::Backend;


		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> 
		decltype(
			FEBackend<Backend>::apply_binary(
							   FEEval<Left, Traits, Backend, 1>::apply(expr.left(), ctx),
				   			   FEEval<Right, Traits, Backend, 1>::apply(expr.right(), ctx),
				   			   expr.operation(),
				   			   ctx
				   			   )
			)
		{
			return FEBackend<Backend>::apply_binary(
				   FEEval<Left, Traits, Backend, 1>::apply(expr.left(), ctx),
	   			   FEEval<Right, Traits, Backend, 1>::apply(expr.right(), ctx),
	   			   expr.operation(),
	   			   ctx
	   			   );
		}
	};

	template<class Left, class Right, class Traits, int IsQuadData>
	class BinaryDelegate<Binary<Left, Right, Minus>, Traits, IsQuadData, 1, 1> {
	public:
		typedef utopia::Binary<Left, Right, Minus> BinExpr;

		static const int Backend = Traits::Backend;
		typedef typename Traits::Vector Vector;

		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> Vector
		{
			return 
				FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx) -
				FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(), ctx);
		}
	};

	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	template<class Left, class Right, class Op, class Traits, int IsQuadData>
	class BinaryDelegate<Binary<Left, Right, Op>, Traits, IsQuadData, 2, 1> {
	public:
		typedef utopia::Binary<Left, Right, Op> BinExpr;

		static const int Backend = Traits::Backend;
		typedef typename Traits::Matrix Matrix;

		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> Matrix
		{
			return FEBackend<Backend>::apply_binary(
				FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx),
				FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(), ctx),
				expr.operation(),
				ctx);
		}
	};


	template<class Left, class Right, class Traits, int IsQuadData>
	class BinaryDelegate<Binary<Left, Right, Plus>, Traits, IsQuadData, 2, 1> {
	public:
		typedef utopia::Binary<Left, Right, Plus> BinExpr;

		static const int Backend = Traits::Backend;
		typedef typename Traits::Matrix Matrix;

		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> Matrix
		{
			return 
				FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx) +
				FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(), ctx);
		}
	};


	template<class Left, class Right, class Traits, int IsQuadData>
	class BinaryDelegate<Binary<Left, Right, Minus>, Traits, IsQuadData, 2, 1> {
	public:
		typedef utopia::Binary<Left, Right, Minus> BinExpr;

		static const int Backend = Traits::Backend;
		typedef typename Traits::Matrix Matrix;

		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> Matrix
		{
			return 
				FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx) -
				FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(), ctx);
		}
	};

	template<class Left, class Right, class Op, class Traits, int Backend, int IsQuadData>
	class FEEval< Binary<Left, Right, Op>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Binary<Left, Right, Op> Expr;
		typedef utopia::BinaryDelegate<Expr, Traits, IsQuadData> DelegateT;

		inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx) -> decltype( DelegateT::apply(expr, ctx) )
		{
			return DelegateT::apply(expr, ctx);
		}
	};


	template<class Left, class Right, class Traits, int Backend, int IsQuadData>
	class FEEval<Binary<Integral<Left>, Right, Plus>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Binary<Integral<Left>, Right, Plus> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	FEEval<Integral<Left>, Traits, Backend, IsQuadData>::apply(expr.left())
	    	)
	    {
	    	return FEEval<Integral<Left>, Traits, Backend, IsQuadData>::apply(expr.left()) + 
	    		   FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right());
	    } 
	};

	template<class Left, class Right, class Traits, int Backend, int IsQuadData>
	class FEEval<Binary<Transposed< Gradient<Left> >, Gradient<Right>, Plus>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Binary<Transposed< Gradient<Left> >, Gradient<Right>, Plus> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEBackend<Backend>::grad_t_plus_grad(1.0, expr.left().expr().expr(), expr.right().expr(), ctx)
	    	    	)
	    {
	    	return FEBackend<Backend>::grad_t_plus_grad(1.0, expr.left().expr().expr(), expr.right().expr(), ctx);
	    } 
	};





	template<class Left, class Right, class Traits, int Backend, int IsQuadData>
	class FEEval<Binary<Gradient<Left>, Transposed< Gradient<Right> >, Plus>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Binary<Gradient<Left>,Transposed< Gradient<Right> >, Plus> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEBackend<Backend>::grad_t_plus_grad(expr.right(), expr.left().expr(), ctx)
	    	    	)
	    {
	    	return FEBackend<Backend>::grad_t_plus_grad(expr.right(), expr.left().expr(), ctx);
	    } 
	};


	template<class Left, class Right, class Traits, int Backend, int IsQuadData>
	class FEEval<Binary<Integral<Left>, Integral<Right>, Plus>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Binary<Integral<Left>, Integral<Right>, Plus> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left().expr(), ctx)
	    	    	)
	    {
	    	return FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left().expr(), ctx) + FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right().expr(), ctx);
	    } 
	};

	template<class Left, class Right, class Traits, int Backend, int IsQuadData>
	class FEEval<Binary<Integral<Left>, Integral<Right>, Minus>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Binary<Integral<Left>, Integral<Right>, Minus> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left().expr(), ctx)
	    	    	)
	    {
	    	return FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left().expr(), ctx) - FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right().expr(), ctx);
	    } 
	};







	template<class T, int Order, class C, class F, class Op, class Traits, int Backend, int IsQuadData>
	class FEEval<Binary<ConstantCoefficient<T, Order>, Interpolate<C, F>, Op>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::ConstantCoefficient<T, Order> Left;
		typedef utopia::Interpolate<C, F> Right;
		typedef utopia::Binary<Left, Right, Op> Expr;

		static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx) -> decltype( FEEval<Right, Traits, Backend, IsQuadData>::apply(expr.right(), ctx) )
		{
			return FEBackend<Backend>::apply_binary(
				expr.left(),
				expr.right(),
				expr.operation(),
				ctx);
		}
	};


	template<class C, class F, class T, int Order, class Op, class Traits, int Backend, int IsQuadData>
	class FEEval<Binary<Interpolate<C, F>, ConstantCoefficient<T, Order>, Op>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Interpolate<C, F> Left;
		typedef utopia::ConstantCoefficient<T, Order> Right;
		typedef utopia::Binary<Left, Right, Op> Expr;

		static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx) -> decltype( FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx) )
		{
			return FEBackend<Backend>::apply_binary(
				expr.left(),
				expr.right(),
				expr.operation(),
				ctx);
		}
	};



	template<class T, int Order, class C, class F, class Traits, int Backend, int IsQuadData>
	class FEEval<Binary<Interpolate<C, F>, ConstantCoefficient<T, Order>, Plus>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::Interpolate<C, F> Left;
		typedef utopia::ConstantCoefficient<T, Order> Right;
		typedef utopia::Binary<Left, Right, Plus> Expr;

		static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx) -> decltype( FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx) )
		{
			return FEBackend<Backend>::apply_binary(
				expr.right(),
				expr.left(),
				expr.operation(),
				ctx);
		}
	};



	template<class Right, class Op, class Traits, int Backend, int IsQuadData>
	class FEEval<Binary<SymbolicTensor<Identity, 2>, Gradient<Right>, Op>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::SymbolicTensor<Identity, 2> Left;
		typedef utopia::Binary<Left, Gradient<Right>, Op> Expr;

		static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx) -> decltype( FEEval<Gradient<Right>, Traits, Backend, IsQuadData>::apply(expr.right(), ctx) )
		{
			return FEBackend<Backend>::apply_binary(
				expr.left(),
				FEEval<Gradient<Right>, Traits, Backend, IsQuadData>::apply(expr.right(), ctx),
				expr.operation(),
				ctx);
		}
	};


	template<class Left, class Op, class Traits, int Backend, int IsQuadData>
	class FEEval<Binary<Left, SymbolicTensor<Identity, 2>, Op>, Traits, Backend, IsQuadData> {
	public:
		typedef utopia::SymbolicTensor<Identity, 2> Right;
		typedef utopia::Binary<Left, Right, Op> Expr;

		static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx) -> decltype( FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx) )
		{
			return FEBackend<Backend>::apply_binary(
				FEEval<Left, Traits, Backend, IsQuadData>::apply(expr.left(), ctx),
				expr.right(),
				expr.operation(),
				ctx);
		}
	};

	////////////////////////////////////////////////////////////////////////////////////////////////


	template<class Left, class Right, class Op, class AssemblyContext>
	class FunctionalTraits<Binary<Left, Right, Op>, AssemblyContext> {
	public:
		inline static int type(const Binary<Left, Right, Op> &expr, AssemblyContext &ctx)  
		{ 
			return std::max(FunctionalTraits<Left, AssemblyContext>::type(expr.left(), ctx), 
							FunctionalTraits<Right, AssemblyContext>::type(expr.right(), ctx));
		}

		inline static int order(const Binary<Left, Right, Op> &expr, const AssemblyContext &ctx) 
		{
			return std::max(FunctionalTraits<Left, AssemblyContext>::order(expr.left(), ctx), 
				   FunctionalTraits<Right, AssemblyContext>::order(expr.right(), ctx));
		}
	};

	template<class Left, class Right, class AssemblyContext>
	class FunctionalTraits<Binary<Left, Right, EMultiplies>, AssemblyContext> {
	public:
		inline static int type(const Binary<Left, Right, EMultiplies> &expr, const AssemblyContext &ctx)  
		{ 
			return std::max(FunctionalTraits<Left, AssemblyContext>::type(expr.left(), ctx), 
							FunctionalTraits<Right, AssemblyContext>::type(expr.right(), ctx));
		}

		inline static int order(const Binary<Left, Right, EMultiplies> &expr, const AssemblyContext &ctx) 
		{
			return FunctionalTraits<Left, AssemblyContext>::order(expr.left(), ctx) +
				   FunctionalTraits<Right, AssemblyContext>::order(expr.right(), ctx);
		}
	};

}


#endif //UTOPIA_FE_EVAL_BINARY_HPP
