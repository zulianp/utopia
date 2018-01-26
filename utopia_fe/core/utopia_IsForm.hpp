#ifndef UTOPIA_IS_FORM_HPP
#define UTOPIA_IS_FORM_HPP

#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_FEIsSubTree.hpp"

namespace utopia {
	///gives false positives
	template<class Expr>
	struct IsForm 
	{	
		typedef utopia::Reduce< Binary<utopia::Any, utopia::Any, EMultiplies>, Plus> Inner;
		typedef Binary<utopia::Any, utopia::Any, Multiplies> Inner2;

		static const int has_inner_product = IsSubTree<Inner, Expr>::value;
		static const int has_scalar_inner_product = IsSubTree<Inner2, Expr>::value;

		static const int has_trial    	   = IsSubTree<TrialFunction<utopia::Any>, Expr>::value;
		static const int has_test          = IsSubTree<TestFunction<utopia::Any>,  Expr>::value;
		static const int has_fun           = has_trial || has_test;
		static const int order             = has_inner_product * (has_trial + has_test);
		static const int value             = has_inner_product && (has_fun);

	};

	template<class Expr>
	struct IsBilinearForm
	{
		static const int value = IsForm<Expr>::order == 2;
	};

	template<class Expr>
	struct IsLinearForm
	{
		static const int value = IsForm<Expr>::order == 2;
	};


	template<class Expr>
	struct IsForm<const Expr> : IsForm<Expr> {};
	
	template<class Expr>
	struct IsForm<const Expr &> : IsForm<Expr> {};

	template<class Expr>
	struct IsForm<Expr &> : IsForm<Expr> {};

	template<class Expr>
	struct IsForm<Expr &&> : IsForm<Expr> {};
}

#endif //UTOPIA_IS_FORM_HPP
