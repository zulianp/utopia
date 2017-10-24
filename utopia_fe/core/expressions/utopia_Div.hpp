#ifndef UTOPIA_FE_DIV_HPP
#define UTOPIA_FE_DIV_HPP

#include "utopia_DifferentialOperator.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

	template<class Expr>
	class Divergence : public DifferentialOperator< Divergence<Expr> > {
	public:
		enum {
			Order = Expr::Order - 1
		};

		typedef typename Expr::Scalar Scalar;

		std::string getClass() const { return "Divergence<" + expr_.getClass() + ">"; }

		inline const Expr &expr() const
		{
			return expr_;
		}

		Divergence(const Expr &expr)
		: expr_(expr)
		{}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
	};

	template<class Derived>
	inline Divergence<Derived> div(const Expression<Derived> &expr) {
		return Divergence<Derived>(expr);
	}

	template<class Expr>
	class Traits< Divergence<Expr> > : public Traits<Expr> {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};
}

#endif //UTOPIA_FE_DIV_HPP
