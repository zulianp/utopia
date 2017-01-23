#ifndef UTOPIA_AUTO_DIFF_EXPR_BINARY_HPP
#define UTOPIA_AUTO_DIFF_EXPR_BINARY_HPP 

#include "utopia_AutoDiffExpr.hpp"

namespace utopia {
	template<class Left, class Right, class Operation>
	class AutoDiffBinary { };

	template<class Left, class Right>
	class AutoDiffBinary<Left, Right, EMultiplies> { 
	public:
		typedef utopia::AutoDiffExpr<Left>  DiffLeft;
		typedef utopia::AutoDiffExpr<Right> DiffRight;

		typedef typename DiffLeft::Type  DLeft;
		typedef typename DiffRight::Type DRight;

		typedef utopia::Binary< Binary<DLeft, Right, EMultiplies>,
								Binary<Left, DRight, EMultiplies>, 
								Plus> ComplexType;

		typedef typename utopia::Simplify<ComplexType>::Type Type;	
								
		inline static UTOPIA_STORE_CONST(Type) make(const Left &left, const Right &right, const EMultiplies &)
		{
			return 	 
					utopia::Simplify<ComplexType>::make
					( 
						e_mul(DiffLeft::make(left), right) +
						e_mul(left, DiffRight::make(right))
					);
		}
	};


	template<class Left, class Right>
	class AutoDiffBinary<Left, Right, Plus> { 
	public:
		typedef utopia::AutoDiffExpr<Left>  DiffLeft;
		typedef utopia::AutoDiffExpr<Right> DiffRight;

		typedef typename DiffLeft::Type  DLeft;
		typedef typename DiffRight::Type DRight;

		typedef utopia::Binary<DLeft, DRight, Plus> ComplexType;
		typedef typename utopia::Simplify<ComplexType>::Type Type;	
								

		inline static UTOPIA_STORE_CONST(Type) make(const Left &left, const Right &right, const Plus &)
		{
			return 	utopia::Simplify<ComplexType>::make
					( 
						DiffLeft::make(left) + DiffRight::make(right)
					);
		}
	};

	template<class Left, class Right, class Operation>
	class AutoDiffExpr< Binary<Left, Right, Operation>, 1> {
	public:

		typedef utopia::AutoDiffBinary<Left, Right, Operation> Diff;
		typedef typename Diff::Type Type;

		inline static UTOPIA_STORE_CONST(Type) make(const Binary<Left, Right, Operation> &expr)
		{
			return Diff::make(expr.left(), expr.right(), expr.operation());
		}
	};


	template<class Left, class Right>
	class AutoDiffBinary<Number<Left>, Right, Multiplies> { 
	public:
		typedef utopia::AutoDiffExpr<Right> DiffRight;
		typedef typename DiffRight::Type DRight;

		typedef utopia::Binary<Number<Left>, DRight, Multiplies> ComplexType;
		typedef typename utopia::Simplify<ComplexType>::Type Type;	

		inline static UTOPIA_STORE_CONST(Type) make(const Left &left, const Right &right, const Multiplies &)
		{
			return 	utopia::Simplify<ComplexType>::make
					( 
						left * DiffRight::make(right)
					);
		}
	};
}

#endif //UTOPIA_AUTO_DIFF_EXPR_BINARY_HPP

