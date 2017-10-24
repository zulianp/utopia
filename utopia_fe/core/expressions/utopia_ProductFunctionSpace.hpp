#ifndef UTOPIA_PRODUCT_FUNCTION_SPACE_HPP
#define UTOPIA_PRODUCT_FUNCTION_SPACE_HPP 

#include "utopia_FormTraits.hpp"
#include "utopia_FunctionSpace.hpp"
#include <utility>

namespace utopia {

	template<class... Spaces>
	class ProductFunctionSpace : public FunctionSpace<ProductFunctionSpace<Spaces...> > {
	public:

		ProductFunctionSpace(const Spaces &...spaces)
		: spaces_(spaces...)
		{ }

		template<int Index>
		inline auto get() -> decltype(std::get<Index>(std::tuple<Spaces...>())) &
		{
			return std::get<Index>(this->spaces_);
		}

	private:
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

	template<class... Spaces>
	inline ProductFunctionSpace<FunctionSpace<Spaces>...> kron_prod(const FunctionSpace<Spaces> &... spaces)
	{
		return ProductFunctionSpace<FunctionSpace<Spaces>...>(spaces...);
	}

	template<class... Spaces>
	class FormTraits< ProductFunctionSpace<Spaces...> > {
	public:
		template<class Space>
		inline static constexpr int order()
		{
			return Traits<Space>::Order;
		}

		template<class First, class...Rest>
		inline static constexpr int order()
		{
			return Traits<First>::Order + order<Rest...>();
		}

		template<class First, class...Rest>
		class GetFirst {
		public:
			typedef First Type;
		};

		typedef typename GetFirst<Spaces...>::Type First;
		typedef typename utopia::FormTraits<First>::Implementation Implementation;
		typedef typename utopia::FormTraits<First>::Scalar Scalar;
		const static int Order = order<Spaces...>();
	};
}

#endif //UTOPIA_PRODUCT_FUNCTION_SPACE_HPP
