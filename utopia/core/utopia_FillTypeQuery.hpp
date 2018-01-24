#ifndef UTOPIA_FILL_TYPE_QUERY_HPP
#define UTOPIA_FILL_TYPE_QUERY_HPP 

#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

	/**********************************************************************************************

	 template<class Expr>
	 class FillTypeQuery {
	 public:
	     enum {
	         value = Expr::FILL_TYPE
	     };
	 };

	 template<class Left, class Right, class Operation>
	 class FillTypeQuery< Binary<Left, Right, Operation> > {
	 public:
	 	enum {
	 			value = ( FillTypeQuery<Left>::value == FillType::DELEGATE ?
	 					  FillTypeQuery<Right>::value  : (
	                       FillTypeQuery<Left>::value == FillType::SCALAR ?
	                            ( FillTypeQuery<Left>::value == FillType::SCALAR ?
	                            	 FillType::SCALAR : FillTypeQuery<Right>::value  ) :
	                    		FillTypeQuery<Left>::value )
	 	              )
	 	};
	 };

	 template<class Left, class Right>
	 class FillTypeQuery< Multiply<Left, Right> > {
	 public:
	 	enum {
	 			value = ( FillTypeQuery<Left>::value == FillType::DELEGATE ?
	 					  FillTypeQuery<Right>::value  : (
	                       FillTypeQuery<Left>::value == FillType::SCALAR ?
	                            ( FillTypeQuery<Left>::value == FillType::SCALAR ?
	                            	 FillType::SCALAR : FillTypeQuery<Right>::value  ) :
	                    		FillTypeQuery<Left>::value )
	 	              )
	 	};
	 };


	 template<class Expr, class Operation>
	 class FillTypeQuery< Unary<Expr, Operation> > {
	 public:
	     enum {
	         value = FillTypeQuery<Expr>::value
	     };
	 };

	 template<class Expr, int Order>
	 class FillTypeQuery< Norm<Expr, Order> > {
	 public:
	     enum {
	         value = FillType::SCALAR
	     };
	 };

	 template<class Expr>
	 class FillTypeQuery< Transposed<Expr> > {
	 public:
	     enum {
	         value = FillTypeQuery<Expr>::value
	     };
	 };

	 template<class Left, class Right>
	 class FillTypeQuery< Construct<Left, Right> > {
	 public:
	     enum {
	         value = FillTypeQuery<Left>::value
	     };
	 };

	 template<class Left, class Right>
	 class FillTypeQuery< Assign<Left, Right> > {
	 public:
	     enum {
	         value = FillTypeQuery<Left>::value
	     };
	 };


	 template<class ResultImpl, class Expr>
	 class SelectFillType {
	 public:
	 	typedef typename utopia::Traits<ResultImpl>::Traits Traits;

	 	enum {
	 		EXPR_FILL_TYPE = FillTypeQuery<Expr>::value
	 	};

	 	enum {
	 		RESULT_FILL_TYPE = Traits::FILL_TYPE
	 	};

	 	enum {
	 		aux_value = (
	 			( RESULT_FILL_TYPE == FillType::SPARSE || ( EXPR_FILL_TYPE == FillType::DELEGATE &&
	 														RESULT_FILL_TYPE != FillType::SCALAR ) )?
	 			RESULT_FILL_TYPE : EXPR_FILL_TYPE )
	 	};

	 	enum {
	 		value = (aux_value == FillType::DELEGATE || aux_value == FillType::SCALAR ) ? FillType::SPARSE : aux_value
	 	};
	 };


	 #define EXPR_AND_FILL_TYPE(Traits, Expr, Result) \
	 	typename utopia::TensorQuery<Traits, Expr::Order, utopia::Traits<Expr>::FILL_TYPE>::Type

    **********************************************************************************************/
}


#endif //UTOPIA_FILL_TYPE_QUERY_HPP
