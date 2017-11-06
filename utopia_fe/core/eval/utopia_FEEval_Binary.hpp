#ifndef UTOPIA_FE_EVAL_BINARY_HPP
#define UTOPIA_FE_EVAL_BINARY_HPP

#include "utopia_FEEval_Empty.hpp"
#include "utopia_FEForwardDeclarations.hpp"

namespace utopia {

	template<class BinExpr, int FormOrder>
	class BinaryDelegate {};

	template<class BinExpr>
	class BinaryDelegate<BinExpr, 0> {
	public:
		// template<class AssemblyContextT>
		// static auto apply(const BinExpr &expr, AssemblyContextT &) -> decltype( Eval<BinExpr>::apply(expr) )
		// {
		// 	return Eval<BinExpr>::apply(expr);
		// }

		template<class AssemblyContextT>
		static auto apply(const BinExpr &expr, AssemblyContextT &) -> const BinExpr &
		{
			return expr;
		}
	};

	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	template<class Left, class Right, class Op>
	class BinaryDelegate<Binary<Left, Right, Op>, 1> {
	public:
		typedef utopia::Binary<Left, Right, Op> BinExpr;

		typedef utopia::Traits<BinExpr> Traits;
		static const int Backend = Traits::Backend;
		typedef typename Traits::Vector Vector;

		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> Vector
		{
			return FEBackend<Backend>::apply_binary(
				FEEval<Left, Traits, Backend>::apply(expr.left(), ctx),
				FEEval<Right, Traits, Backend>::apply(expr.right(), ctx),
				ctx);
		}
	};


	template<class Left, class Right>
	class BinaryDelegate<Binary<Left, Right, Plus>, 1> {
	public:
		typedef utopia::Binary<Left, Right, Plus> BinExpr;

		typedef utopia::Traits<BinExpr> Traits;
		static const int Backend = Traits::Backend;
		typedef typename Traits::Vector Vector;

		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> Vector
		{
			return FEBackend<Backend>::apply_binary(
				FEEval<Left, Traits, Backend>::apply(expr.left(), ctx),
				FEEval<Right, Traits, Backend>::apply(expr.right(), ctx),
				ctx);
		}
	};

	template<class Left, class Right>
	class BinaryDelegate<Binary<Left, Right, Minus>, 1> {
	public:
		typedef utopia::Binary<Left, Right, Minus> BinExpr;

		typedef utopia::Traits<BinExpr> Traits;
		static const int Backend = Traits::Backend;
		typedef typename Traits::Vector Vector;

		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> Vector
		{
			return FEBackend<Backend>::apply_binary(
				FEEval<Left, Traits, Backend>::apply(expr.left(), ctx),
				FEEval<Right, Traits, Backend>::apply(expr.right(), ctx),
				ctx);
		}
	};

	//////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////

	template<class Left, class Right, class Op>
	class BinaryDelegate<Binary<Left, Right, Op>, 2> {
	public:
		typedef utopia::Binary<Left, Right, Op> BinExpr;

		typedef utopia::Traits<BinExpr> Traits;
		static const int Backend = Traits::Backend;
		typedef typename Traits::Matrix Matrix;

		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> Matrix
		{
			return FEBackend<Backend>::apply_binary(
				FEEval<Left, Traits, Backend>::apply(expr.left(), ctx),
				FEEval<Right, Traits, Backend>::apply(expr.right(), ctx),
				ctx);
		}
	};


	template<class Left, class Right>
	class BinaryDelegate<Binary<Left, Right, Plus>, 2> {
	public:
		typedef utopia::Binary<Left, Right, Plus> BinExpr;

		typedef utopia::Traits<BinExpr> Traits;
		static const int Backend = Traits::Backend;
		typedef typename Traits::Matrix Matrix;

		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> Matrix
		{
			return 
				FEEval<Left, Traits, Backend>::apply(expr.left(), ctx) +
				FEEval<Right, Traits, Backend>::apply(expr.right(), ctx);
		}
	};


	template<class Left, class Right>
	class BinaryDelegate<Binary<Left, Right, Minus>, 2> {
	public:
		typedef utopia::Binary<Left, Right, Minus> BinExpr;

		typedef utopia::Traits<BinExpr> Traits;
		static const int Backend = Traits::Backend;
		typedef typename Traits::Matrix Matrix;

		static auto apply(const BinExpr &expr, AssemblyContext<Backend> &ctx) -> Matrix
		{
			return 
				FEEval<Left, Traits, Backend>::apply(expr.left(), ctx) +
				FEEval<Right, Traits, Backend>::apply(expr.right(), ctx);
		}
	};

	template<class Left, class Right, class Op, class Traits, int Backend>
	class FEEval< Binary<Left, Right, Op>, Traits, Backend> {
	public:
		typedef utopia::Binary<Left, Right, Op> Expr;
		static const int has_trial = IsSubTree<TrialFunction<utopia::Any>, Expr>::value;
		static const int has_test  = IsSubTree<TestFunction<utopia::Any>,  Expr>::value;
		typedef utopia::BinaryDelegate<Expr, has_trial + has_test> DelegateT;

		inline static auto apply(const Expr &expr, AssemblyContext<Backend> &ctx) -> decltype( DelegateT::apply(expr, ctx) )
		{
			return DelegateT::apply(expr, ctx);
		}
	};


	template<class Left, class Right, class Traits, int Backend>
	class FEEval<Binary<Integral<Left>, Right, Plus>, Traits, Backend> {
	public:
		typedef utopia::Binary<Integral<Left>, Right, Plus> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	FEEval<Integral<Left>, Traits, Backend>::apply(expr.left())
	    	)
	    {
	    	return FEEval<Integral<Left>, Traits, Backend>::apply(expr.left()) + 
	    		   FEEval<Right, Traits, Backend>::apply(expr.right());
	    } 
	};

	template<class Left, class Space, class Op, class Traits, int Backend>
	class FEEval<Binary<Left, TrialFunction<Space>, Op>, Traits, Backend> {
	public:
		typedef Binary<Left, TrialFunction<Space>, Op> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    		FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), expr.operation(), ctx)
	    	)
	    {
	    	return FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), expr.operation(), ctx);
	    } 
	};

	template<class Left, class Space, class Op, class Traits, int Backend>
	class FEEval<Binary<Left, TestFunction<Space>, Op>, Traits, Backend> {
	public:
		typedef Binary<Left, TestFunction<Space>, Op> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), expr.operation(), ctx)
	    	    	)
	    {
	    	return FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), expr.operation(), ctx);
	    } 
	};

	template<class Left, class Function, class Op, class Traits, int Backend>
	class FEEval<Binary<Left, Gradient<Function>, Op>, Traits, Backend> {
	public:
		typedef Binary<Left, Gradient<Function>, Op> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), expr.operation(), ctx)
	    	    	)
	    {
	    	return FEBackend<Backend>::apply_binary(FEEval<Left, Traits, Backend>::apply(expr.left(), ctx), expr.right(), expr.operation(), ctx);
	    } 
	};

	template<class Left, class Right, class Traits, int Backend>
	class FEEval<Binary<Transposed< Gradient<Left> >, Gradient<Right>, Plus>, Traits, Backend> {
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

	template<typename T, class Left, class Right, class Traits, int Backend>
	class FEEval< Binary<Number<T>, 
						 Binary<Transposed< Gradient<Left> >, Gradient<Right>, Plus>,
						 Multiplies>,
				  Traits, Backend> {
	public:
		typedef utopia::Binary<Number<T>, 
						 Binary<Transposed< Gradient<Left> >, Gradient<Right>, Plus>,
						 Multiplies> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEBackend<Backend>::grad_t_plus_grad(expr.left(), expr.right().left().expr().expr(), expr.right().right().expr(), ctx)
	    	    	)
	    {
	    	return FEBackend<Backend>::grad_t_plus_grad(expr.left(), expr.right().left().expr().expr(), expr.right().right().expr(), ctx);
	    } 
	};



	template<class Left, class Right, class Traits, int Backend>
	class FEEval<Binary<Gradient<Left>, Transposed< Gradient<Right> >, Plus>, Traits, Backend> {
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


	template<class Left, class Right, class Traits, int Backend>
	class FEEval<Binary<Integral<Left>, Integral<Right>, Plus>, Traits, Backend> {
	public:
		typedef utopia::Binary<Integral<Left>, Integral<Right>, Plus> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEEval<Left, Traits, Backend>::apply(expr.left().expr(), ctx)
	    	    	)
	    {
	    	return FEEval<Left, Traits, Backend>::apply(expr.left().expr(), ctx) + FEEval<Right, Traits, Backend>::apply(expr.right().expr(), ctx);
	    } 
	};

	template<class Left, class Right, class Traits, int Backend>
	class FEEval<Binary<Integral<Left>, Integral<Right>, Minus>, Traits, Backend> {
	public:
		typedef utopia::Binary<Integral<Left>, Integral<Right>, Minus> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEEval<Left, Traits, Backend>::apply(expr.left().expr(), ctx)
	    	    	)
	    {
	    	return FEEval<Left, Traits, Backend>::apply(expr.left().expr(), ctx) - FEEval<Right, Traits, Backend>::apply(expr.right().expr(), ctx);
	    } 
	};


	template<class Left, class Right, class Traits, int Backend>
	class FEEval<Binary<Number<Left>, Integral<Right>, Multiplies>, Traits, Backend> {
	public:
		typedef utopia::Binary<Number<Left>, Integral<Right>, Multiplies> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEEval<Right, Traits, Backend>::apply(expr.right().expr(), ctx)
	    	    	)
	    {
	    	return static_cast<Left>(expr.left()) * FEEval<Right, Traits, Backend>::apply(expr.right().expr(), ctx);
	    } 
	};


	template<class Left, class Right, class Op, class Traits, int Backend>
	class FEEval<Binary<Number<Left>, Reduce<Right, Op>, Multiplies>, Traits, Backend> {
	public:
		typedef utopia::Binary<Number<Left>, Reduce<Right, Op>, Multiplies> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEEval<Reduce<Right, Op>, Traits, Backend>::apply(expr.right(), ctx)
	    	    	)
	    {
	    	return static_cast<Left>(expr.left()) * FEEval<Reduce<Right, Op>, Traits, Backend>::apply(expr.right(), ctx);
	    } 
	};

	template<class Left, class Right, class Traits, int Backend>
	class FEEval<Binary<Integral<Left>, Number<Right>, Multiplies>, Traits, Backend> {
	public:
		typedef utopia::Binary<Integral<Left>, Number<Right>, Multiplies> Expr;


	    inline static auto apply(
	    	const Expr &expr,
	    	AssemblyContext<Backend> &ctx) -> decltype( 
	    	    		FEEval<Left, Traits, Backend>::apply(expr.left(), ctx)
	    	    	)
	    {
	    	return FEEval<Left, Traits, Backend>::apply(expr.left(), ctx) * FEEval<Right, Traits, Backend>::apply(expr.right(), ctx);
	    } 
	};
	////////////////////////////////////////////////////////////////////////////////////////////////


	template<class Left, class Right, class Op, class AssemblyContext>
	class FunctionalTraits<Binary<Left, Right, Op>, AssemblyContext> {
	public:
		inline static int type(const Binary<Left, Right, Op> &expr, const AssemblyContext &ctx)  
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


	template<class Left, class Right, class AssemblyContext>
	class FunctionalTraits<Binary<Left, Right, Multiplies>, AssemblyContext> {
	public:
		inline static int type(const Binary<Left, Right, Multiplies> &expr, const AssemblyContext &ctx)  
		{ 
			return std::max(FunctionalTraits<Left, AssemblyContext>::type(expr.left(), ctx), 
							FunctionalTraits<Right, AssemblyContext>::type(expr.right(), ctx));
		}

		inline static int order(const Binary<Left, Right, Multiplies> &expr, const AssemblyContext &ctx) 
		{
			return FunctionalTraits<Left, AssemblyContext>::order(expr.left(), ctx) +
				   FunctionalTraits<Right, AssemblyContext>::order(expr.right(), ctx);
		}
	};
}


#endif //UTOPIA_FE_EVAL_BINARY_HPP
