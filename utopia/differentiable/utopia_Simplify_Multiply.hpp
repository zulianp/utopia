#ifndef DO_NOT_SIMPLIFY_DERIVATIVES

#ifndef UTOPIA_SIMPLIFY_MULTIPLY_HPP
#define UTOPIA_SIMPLIFY_MULTIPLY_HPP

namespace utopia {

	template<class Left, class Right>
	class MultiplySimplify {
	public:
		typedef utopia::Multiply<Left, Right> Type;

		inline static const Type &make(const Type &expr)
		{
			return expr;
		}
	};

	template<class Expr, int Order>
	class MultiplySimplify<Expr, Factory<Identity, Order> > {
	public:
		typedef typename Simplify<Expr>::Type Type;
		
		static UTOPIA_STORE_CONST(Type) make(const Multiply<Expr, Factory<Identity, Order> > &expr)
		{
			return Simplify<Expr>::make(expr.left());
		}
	};
	
	template<class Expr, int Order>
	class MultiplySimplify< Factory<Identity, Order>, Expr > {
	public:
		typedef typename Simplify<Expr>::Type Type;
		
		static UTOPIA_STORE_CONST(Type) make(const Multiply<Factory<Identity, Order>, Expr>  &expr)
		{
			return Simplify<Expr>::make(expr.right());
		}
	};
	
	template<class Expr, int Order>
	class MultiplySimplify< Factory<Zeros, Order>, Expr >{
	public:
		typedef utopia::Factory<Zeros, Order> Type;
		
		static UTOPIA_STORE_CONST(Type) make(const Multiply< Factory<Zeros, Order>, Expr>  &expr)
		{
			return expr.left();
		}
	};
	
	template<class Expr, int Order>
	class MultiplySimplify<Expr, Factory<Zeros, Order> > {
	public:
		typedef utopia::Factory<Zeros, Order> Type;
		
		static UTOPIA_STORE_CONST(Type) make(const Multiply< Expr, Factory<Zeros, Order> >  &expr)
		{
			return expr.right();
		}
	};

	template<class Left, class Right>
	class Simplify< Multiply<Left, Right> > {
	public:
		typedef typename utopia::Simplify<Left>::Type  SLeft;
		typedef typename utopia::Simplify<Right>::Type SRight;
		
		typedef utopia::MultiplySimplify<SLeft, SRight> Sim;
		typedef typename Sim::Type Type;

		static UTOPIA_STORE_CONST(Type) make(const Multiply<Left, Right> &expr)
		{
			return Sim::make(
						utopia::Simplify<Left>::make(expr.left()) * 
						utopia::Simplify<Right>::make(expr.right()) );
		}
	};
}

#endif //UTOPIA_SIMPLIFY_MULTIPLY_HPP
#endif //DO_NOT_SIMPLIFY_DERIVATIVES
