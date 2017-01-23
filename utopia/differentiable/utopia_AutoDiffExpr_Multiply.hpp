#ifndef UTOPIA_AUTO_DIFF_EXPR_MULTIPLY_HPP
#define UTOPIA_AUTO_DIFF_EXPR_MULTIPLY_HPP

#include "utopia_Differentiable.hpp"
#include "utopia_Simplify.hpp"
#include "utopia_AutoDiffExpr.hpp"

namespace utopia {
	
	template<class Left, class Right>
	class AutoDiffExpr< Multiply<Left, Right>, 1> {
	public:
		typedef typename utopia::AutoDiffExpr<Left>  DiffLeft;
		typedef typename utopia::AutoDiffExpr<Right> DiffRight;
		
		typedef typename DiffLeft::Type  DLeft;
		typedef typename DiffRight::Type DRight;
		
		typedef utopia::Binary< utopia::Multiply<DLeft, Right>, 
							    utopia::Multiply<Left, DRight>, Plus > ComplexType;
		
		typedef typename utopia::Simplify<ComplexType>::Type Type;
		
		static UTOPIA_STORE_CONST(Type) make(const Multiply<Left, Right> &expr)
		{
			return utopia::Simplify<ComplexType>::make(
													   		DiffLeft::make(expr.left()) * expr.right() +
													   		expr.left() * DiffRight::make(expr.right())
													   );
		}
	};	
}

#endif //UTOPIA_AUTO_DIFF_EXPR_MULTIPLY_HPP
