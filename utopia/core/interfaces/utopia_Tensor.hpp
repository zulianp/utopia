#ifndef UTOPIA_TENSOR_HPP
#define UTOPIA_TENSOR_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Expression.hpp"
#include "utopia_Operators.hpp"
#include "utopia_Select.hpp"
// #include "utopia_Eval.hpp"

namespace utopia {

	//CRTP type
	template<typename Derived, int Order_>
	class Tensor : 
	// public Expression<Derived> {
	public Expression<Tensor<Derived, Order_>>,
	public Selectable<Tensor<Derived, Order_>> {
	public:
		static const int Order   = Order_;
		static const int StoreAs = UTOPIA_BY_REFERENCE;

		using Scalar = typename Traits<Derived>::Scalar;

		CONST_DERIVED_CRT(Derived);
		DERIVED_CRT(Derived);

		Tensor() {}

		virtual void construct(const Derived &other) { assign(other); }
		virtual void construct(Derived &&other) { assign(other); }

		virtual void assign(const Derived &other) = 0;
		virtual void assign(Derived &&other) = 0;

		template<class Expr>
		void construct_eval(const Expression<Expr> &expr)
		{
			using C = utopia::Construct<Tensor, Expr>;
			Eval<C, Traits<Derived>, Traits<Derived>::Backend>::apply(C(*this, expr.derived()));
		}

		template<class Expr>
		void assign_eval(const Expression<Expr> &expr)
		{
			using A = utopia::Assign<Tensor, Expr>;
			Eval<A, Traits<Derived>, Traits<Derived>::Backend>::apply(A(*this, expr.derived()));
		}

		// template<class Expr>
		void assign_eval(const Tensor &expr)
		{
			assign(expr.derived());
		}

		Tensor &operator=(const Tensor &t) 
		{
			if(this == &t) return *this;
			assign(t);
			return *this;
		}

		Tensor &operator=(Tensor &&t) 
		{
			if(this == &t) return *this;
			assign(t);
			return *this;
		}

		template<class Expr>
		Derived &operator*=(const Expression<Expr> &expr)
		{
			using I = utopia::InPlace<Tensor, Expr, Multiplies>;
			Eval<I, Traits<Derived>, Traits<Derived>::Backend>::apply(I(*this, expr.derived()));
		    return derived();
		}

		Derived &operator*=(const Scalar value)
		{
			using I = utopia::InPlace<Tensor, Number<Scalar>, Multiplies>;
			Eval<I, Traits<Derived>, Traits<Derived>::Backend>::apply(I(*this, value));
		    return derived();
		}

		template<class Expr>
		Derived &operator+=(const Expression<Expr> &expr)
		{
			using I = utopia::InPlace<Tensor, Expr, Plus>;
			Eval<I, Traits<Derived>, Traits<Derived>::Backend>::apply(I(*this, expr.derived()));
		    return derived();
		}

		template<class Expr>
		Derived &operator-=(const Expression<Expr> &expr)
		{
			using I = utopia::InPlace<Tensor, Expr, Minus>;
			Eval<I, Traits<Derived>, Traits<Derived>::Backend>::apply(I(*this, expr.derived()));
		    return derived();
		}

		template<class Expr>
		Derived &operator/=(const Expression<Expr> &expr)
		{
			using I = utopia::InPlace<Tensor, Expr, Divides>;
			Eval<I, Traits<Derived>, Traits<Derived>::Backend>::apply(I(*this, expr.derived()));
		    return derived();
		}

	};

	template<typename Derived, int Order_>
	class Traits<Tensor<Derived, Order_>> : public Traits<Derived> {
	public:
	    static const int Order = Order;
	};

	// template<class Derived>
	// inline Size size(const Tensor<Derived, 2> &t)
	// {
	//     return t.derived().size();
	// }

	// template<class Derived>
	// inline Size local_size(const Tensor<Derived, 2> &t)
	// {
	//     return t.derived().local_size();
	// }

	template<class Derived, int Order>
	inline auto size(const Tensor<Derived, Order> &t)
	{
	    return t.derived().size();
	}

	template<class Derived, int Order>
	inline auto local_size(const Tensor<Derived, Order> &t)
	{
	    return t.derived().local_size();
	}

}

#endif //UTOPIA_TENSOR_HPP
