#ifndef UTOPIA_GRAD_HPP
#define UTOPIA_GRAD_HPP 

#include "utopia_DifferentialOperator.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"
#include "utopia_FunctionalTraits.hpp"

namespace utopia {
	template<class Expr>
	class Gradient : public DifferentialOperator< Gradient<Expr> > {
	public:
		static const int Order = Expr::Order + 1;
		typedef typename Expr::Scalar Scalar;

		std::string getClass() const override { return "Gradient<" + expr_.getClass() + ">"; }

		inline const Expr &expr() const
		{
			return expr_;
		}

		Gradient(const Expr &expr)
		: expr_(expr)
		{}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
	};

	template<class Derived>
	inline Gradient<Derived> grad(const Expression<Derived> &expr) {
		return Gradient<Derived>(expr.derived());
	}

	template<class Expr>
	class Traits< Gradient<Expr> > : public Traits<Expr> {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};

	template<class Expr, class AssemblyContext>
	class FunctionalTraits<Gradient<Expr>, AssemblyContext>  {
	public:
		inline static int type(const Gradient<Expr> &expr,  const AssemblyContext &ctx) { return FunctionalTraits<Expr, AssemblyContext>::type(expr.expr(), ctx);  }
		inline static int order(const Gradient<Expr> &expr, const AssemblyContext &ctx) { return FunctionalTraits<Expr, AssemblyContext>::order(expr.expr(), ctx) - 1; }
	};

}

#endif //UTOPIA_GRAD_HPP
