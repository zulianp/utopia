#ifndef UTOPIA_TIME_DERIVATIVE_HPP
#define UTOPIA_TIME_DERIVATIVE_HPP 

#include "utopia_StoreAs.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
	template<class Expr_, int OrderOfDifferentiation_>
	class TimeDerivative : public Expression< TimeDerivative<Expr_, OrderOfDifferentiation_> > {
	public:

		typedef Expr_ Expr;
		
		enum {
			Order = Expr::Order
		};

		enum {
			OrderOfDifferentiation = OrderOfDifferentiation_
		};

		typedef typename Expr::Scalar Scalar;

		TimeDerivative(const Expr &expr)
		: expr_(expr)
		{}

		std::string getClass() const { return "TimeDerivative<" + expr_.getClass() + ">"; }

		
		inline const Expr &expr() const
		{
			return expr_;
		}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
	};

	template<class Derived>
	TimeDerivative<Derived, 1> dt(const Expression<Derived> &expr)
	{
		return expr.derived();
	}

	template<class Expr, int OrderOfDifferentiation>
	class Traits< TimeDerivative<Expr, OrderOfDifferentiation> > : public Traits<Expr> {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};

	template<typename T>
	class DeltaT : public Number<T> {
	public:
		DeltaT(const T dt = 0.1)
		: Number<T>(dt){}
	};

	template<typename T>
	class Traits< DeltaT<T> > {
	public:
		enum {
			FILL_TYPE = utopia::FillType::SCALAR
		};
	};

}

#endif //UTOPIA_TIME_DERIVATIVE_HPP
