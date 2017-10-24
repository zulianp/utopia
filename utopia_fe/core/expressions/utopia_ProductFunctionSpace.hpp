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
}

#endif //UTOPIA_PRODUCT_FUNCTION_SPACE_HPP
