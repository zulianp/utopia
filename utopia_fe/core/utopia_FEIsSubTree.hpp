#ifndef UTOPIA_FE_IS_SUBTREE_HPP
#define UTOPIA_FE_IS_SUBTREE_HPP

#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_IsSubTree.hpp"
#include "utopia_Any.hpp"

namespace utopia {

	//trial function
	template<class Space>
	class IsSubTree< TrialFunction<Space>, TrialFunction<Space> > {
	public:
		static const int value = 1;
	};

	template<class Expr, class Space>
	class IsSubTree< Expr, TrialFunction<Space> > {
	public:
		static const int value = 0;
	};

	template<class Space>
	class IsSubTree< TrialFunction<utopia::Any>, TrialFunction<Space> > {
	public:
		static const int value = 1;
	};

	//test function
	template<class Space>
	class IsSubTree< TestFunction<Space>, TestFunction<Space> > {
	public:
		static const int value = 1;
	};

	template<class Expr, class Space>
	class IsSubTree< Expr, TestFunction<Space> > {
	public:
		static const int value = 0;
	};

	template<class Space>
	class IsSubTree< TestFunction<utopia::Any>, TestFunction<Space> > {
	public:
		static const int value = 1;
	};

	//gradient
	template<class Inner>
	class IsSubTree< Gradient<Inner>, Gradient<Inner> > {
	public:
		static const int value = 1;
	};

	template<class Expr, class Inner>
	class IsSubTree< Expr, Gradient<Inner> > {
	public:
		static const int value = IsSubTree<Expr, Inner>::value;
	};

	//divergence
	template<class Inner>
	class IsSubTree< Divergence<Inner>, Divergence<Inner> > {
	public:
		static const int value = 1;
	};

	template<class Expr, class Inner>
	class IsSubTree< Expr, Divergence<Inner> > {
	public:
		static const int value = IsSubTree<Expr, Inner>::value;
	};

	//curl
	template<class Inner>
	class IsSubTree< Curl<Inner>, Curl<Inner> > {
	public:
		static const int value = 1;
	};

	template<class Expr, class Inner>
	class IsSubTree< Expr, Curl<Inner> > {
	public:
		static const int value = IsSubTree<Expr, Inner>::value;
	};

	//time-derivative
	template<class Inner, int Order>
	class IsSubTree< TimeDerivative<Inner, Order>, TimeDerivative<Inner, Order> > {
	public:
		static const int value = 1;
	};

	template<class Expr, class Inner, int Order>
	class IsSubTree< Expr, TimeDerivative<Inner, Order> > {
	public:
		static const int value = IsSubTree<Expr, Inner>::value;
	};

	//integral
	template<class Inner>
	class IsSubTree< Integral<Inner>, Integral<Inner> > {
	public:
		static const int value = 1;
	};

	template<class Expr, class Inner>
	class IsSubTree< Expr, Integral<Inner> > {
	public:
		static const int value = IsSubTree<Expr, Inner>::value;
	};


	template<class Inner>
	class IsSubTree< Integral<utopia::Any>, Integral<Inner> > {
	public:
		static const int value = 1;
	};

	///rest of specializations
}

#endif //UTOPIA_FE_IS_SUBTREE_HPP
