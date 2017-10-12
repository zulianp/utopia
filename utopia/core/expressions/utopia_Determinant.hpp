#ifndef UTOPIA_DETERMINANT_HPP
#define UTOPIA_DETERMINANT_HPP 

#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

	template<class Expr>
	class Determinant : public Expression< Determinant<Expr> > {
	public:
		static_assert(Expr::Order >= 2, "must be a 2nd order tensor or greater");

		static const int Order = 0;

		inline const Expr &expr() const 
		{
			return expr_;
		}

		operator typename Traits<Determinant>::Scalar() const
		{
		    Evaluator<typename Traits<Determinant>::Vector, Traits<Determinant>::Backend> eval;
		    return eval.eval(*this);
		}

		Determinant(const Expr &expr) : expr_(expr) {}

	private:
		UTOPIA_STORE_CONST(Expr) expr_;
	};


	template<class Expr>
	class Traits< Determinant<Expr> > : public Traits<Expr> {};

	template<class Expr>
	inline Size size(const Determinant<Expr> & /*expr*/)
	{
		return {1};
	}

	 /**
     * @ingroup     reductions
     * @brief       \f$ tr(A) = \sum_{i = 0}^{n = dim(A)} A_{ii} \f$. \n
     *
     */
	template<class Derived>
	Determinant<Derived> det(const Expression<Derived> &expr)
	{
	    static_assert(Derived::Order == 2, "expr must be a tensor of order 2");
	    return expr.derived();
	}
}

#endif  //UTOPIA_DETERMINANT_HPP
