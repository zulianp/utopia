#ifndef UTOPIA_FE_CURL_HPP
#define UTOPIA_FE_CURL_HPP 

#include "utopia_DifferentialOperator.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
	template<class Expr>
	class Curl : public DifferentialOperator< Curl<Expr> > {
	public:
		enum {
			Order = Expr::Order - 1
		};

		typedef typename Expr::Scalar Scalar;

		std::string getClass() const override { return "Curl<" + expr_.getClass() + ">"; }

		inline const Expr &expr() const
		{
			return expr_;
		}

		Curl(const Expr &expr)
		: expr_(expr)
		{}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
	};

	template<class Derived>
	inline Curl<Derived> curl(const Expression<Derived> &expr) {
		return Curl<Derived>(expr);
	}

	template<class Expr>
	class Traits< Curl<Expr> > : public Traits<Expr> {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};
}

#endif //UTOPIA_FE_CURL_HPP
