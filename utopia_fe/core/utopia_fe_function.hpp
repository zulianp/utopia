#ifndef UTOPIA_FE_FUNCTION_HPP
#define UTOPIA_FE_FUNCTION_HPP 

#include "utopia_Expression.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Base.hpp"
#include <utility>

namespace utopia {

	static const int EMPTY_BACKEND_FLAG = -3;

	class Empty {};

	template<>
	class Traits<Empty> {
	public:
		static const int Backend = EMPTY_BACKEND_FLAG;
		static const int Order = 1;
		typedef double Scalar;

	};

	template<class T>
	class FormTraits {};

	template<class Derived>
	class FunctionSpace {
	public:
		DERIVED_CRT(Derived);
		CONST_DERIVED_CRT(Derived);
	};


	template<class... Spaces>
	class ProductFunctionSpace : public FunctionSpace<ProductFunctionSpace<Spaces...> > {
	public:

		ProductFunctionSpace(const Spaces &...spaces)
		: spaces_(spaces...)
		{ }

	private:
		// UTOPIA_STORE(LeftSpace) left_;
		// UTOPIA_STORE(RightSpace) right_;
		std::tuple<Spaces...> spaces_;
	};

	template<class Left, class Right>
	class FormTraits<ProductFunctionSpace<Left, Right>> {
	public:
		static const int Backend = FormTraits<Left>::Backend;
		static const int Order   = FormTraits<Left>::Order + FormTraits<Right>::Order;
		typedef typename utopia::FormTraits<Left>::Implementation Implementation;
		typedef typename utopia::FormTraits<Left>::Scalar Scalar;
	};

	template<class LeftSpace, class RightSpace>
	inline ProductFunctionSpace<LeftSpace, RightSpace> operator * (const FunctionSpace<LeftSpace> &left, const FunctionSpace<RightSpace> &right)
	{
		return ProductFunctionSpace<LeftSpace, RightSpace>(left.derived(), right.derived());
	}

	template<class Derived>
	class BasisFunction : public Expression<Derived> {};


	template<class FunctionSpaceT>
	class TrialFunction : public BasisFunction< TrialFunction<FunctionSpaceT> > {
	public:
		static const int Order = utopia::FormTraits<FunctionSpaceT>::Order;
		typedef typename utopia::FormTraits<FunctionSpaceT>::Scalar Scalar;
		typedef typename utopia::FormTraits<FunctionSpaceT>::Implementation Implementation;
	};

	template<class FunctionSpaceT>
	class TestFunction : public BasisFunction< TestFunction<FunctionSpaceT> > {
	public:
		static const int Order = utopia::FormTraits<FunctionSpaceT>::Order;
		typedef typename utopia::FormTraits<FunctionSpaceT>::Scalar Scalar;
		typedef typename utopia::FormTraits<FunctionSpaceT>::Implementation Implementation;
	};

	class EmptyFunctionSpace : public FunctionSpace<EmptyFunctionSpace> {
	public:
	};

	template<>
	class FormTraits<EmptyFunctionSpace> : public Traits<Empty> {
	public:
		typedef utopia::Empty Implementation;
	};

	template<class Derived>
	TrialFunction<Derived> trial(const FunctionSpace<Derived> &space)
	{
		//FIXME
		return TrialFunction<Derived>();
	}

	template<class Derived>
	TestFunction<Derived> test(const FunctionSpace<Derived> &space)
	{
		//FIXME
		return TestFunction<Derived>();
	}
}

#endif //UTOPIA_FE_FUNCTION_HPP
