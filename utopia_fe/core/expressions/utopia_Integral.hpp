#ifndef UTOPIA_INTEGRAL_HPP
#define UTOPIA_INTEGRAL_HPP 

#include "utopia_Expression.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_FEExpression.hpp"
#include "utopia_FEIsSubTree.hpp"

namespace utopia {

	template<class Expr_>
	class Integral : public Expression< Integral<Expr_> >, public FEExpression {
	public:
		typedef Expr_ Expr;
		static const int Order = Expr::Order;


		typedef typename Expr::Scalar Scalar;

		std::string getClass() const override { return "Integral<" + expr_.getClass() + ">"; }

		Integral(const Expr &expr, const int block_id = -1)
		: expr_(expr), block_id_(block_id)
		{}

		inline const Expr &expr() const
		{
			return expr_;
		}

		inline int block_id() const
		{
			return block_id_;
		}

		inline int has_block_id() const
		{
			return block_id_ != -1;
		}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
		int block_id_;
	};


	template<class Derived>
	inline Integral<Derived> integral(const Expression<Derived> &expr) {
		static_assert(!IsSubTree<Integral<utopia::Any>, Derived>::value, "nested integrals are not allowed");
		return Integral<Derived>(expr);
	}

	template<class Derived>
	inline Integral<Derived> integral(const Expression<Derived> &expr, const int block_id) {
		static_assert(!IsSubTree<Integral<utopia::Any>, Derived>::value, "nested integrals are not allowed");
		return Integral<Derived>(expr, block_id);
	}

	template<class Expr>
	class Traits< Integral<Expr> > : public Traits<Expr> {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};

	template<class Expr, class AssemblyContext>
	class FunctionalTraits<Integral<Expr>, AssemblyContext>  {
	public:
		inline static int type(const Integral<Expr> &expr,  const AssemblyContext &ctx) { return FunctionalTraits<Expr, AssemblyContext>::type(expr.expr(), ctx);  }
		inline static int order(const Integral<Expr> &expr, const AssemblyContext &ctx) { return FunctionalTraits<Expr, AssemblyContext>::order(expr.expr(), ctx); }
	};
}

#endif //UTOPIA_INTEGRAL_HPP
