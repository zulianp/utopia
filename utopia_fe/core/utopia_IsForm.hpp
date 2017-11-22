#ifndef UTOPIA_IS_FORM_HPP
#define UTOPIA_IS_FORM_HPP

#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_FEIsSubTree.hpp"

namespace utopia {
	template<class Expr>
	struct IsBilinearForm
	{
		static const int has_trial = IsSubTree<TrialFunction<utopia::Any>, Expr>::value;
		static const int has_test  = IsSubTree<TestFunction<utopia::Any>,  Expr>::value;
		static const int value = (has_trial + has_test) == 2;
	};

	template<class Expr>
	struct IsLinearForm
	{
		static const int has_trial = IsSubTree<TrialFunction<utopia::Any>, Expr>::value;
		static const int has_test  = IsSubTree<TestFunction<utopia::Any>,  Expr>::value;
		static const int value = (has_trial + has_test) == 1;
	};
}

#endif //UTOPIA_IS_FORM_HPP
