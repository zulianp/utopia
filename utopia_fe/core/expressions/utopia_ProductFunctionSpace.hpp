#ifndef UTOPIA_PRODUCT_FUNCTION_SPACE_HPP
#define UTOPIA_PRODUCT_FUNCTION_SPACE_HPP 

#include "utopia_FormTraits.hpp"
#include "utopia_FunctionSpace.hpp"
#include <utility>
#include <tuple>

namespace utopia {

	template<class Space, int Begin, int End>
	class SubspaceIterator {
	public:
		template<class Fun>
		void visit(Fun fun)
		{
			fun(Begin, space.template get<Begin>() );
			SubspaceIterator<Space, Begin + 1, End> next(space);
			next.visit(fun);
		}

		SubspaceIterator(Space &space) : space(space) {}
		Space &space;
	};

	template<class Space, int Begin>
	class SubspaceIterator<Space, Begin, Begin + 1> {
	public:
		template<class Fun>
		void visit(Fun fun)
		{
			fun(Begin, space.template get<Begin>());
		}

		SubspaceIterator(Space &space) : space(space) {}
		Space &space;
	};

	template<class... Spaces>
	class ProductFunctionSpace : public FunctionSpace<ProductFunctionSpace<Spaces...> > {
	public:

		static const int NumberOfSubSpaces = std::tuple_size< std::tuple<Spaces...> >::value;

		ProductFunctionSpace(const Spaces &...spaces)
		: spaces_(spaces...)
		{ }

		template<int Index>
		inline auto get() -> decltype(std::get<Index>(std::tuple<Spaces...>())) &
		{
			return std::get<Index>(this->spaces_);
		}

		inline constexpr static std::size_t size() 
		{
			return NumberOfSubSpaces;
		}

		template<class Fun>
		void each(Fun fun)
		{
			SubspaceIterator<ProductFunctionSpace, 0, NumberOfSubSpaces> iter(*this);
			iter.visit(fun);
		}

	private:
		std::tuple<Spaces...> spaces_;
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

	//base case
	template<class Left, class Right>
	class Traits<ProductFunctionSpace<Left, Right>> {
	public:
		typedef typename utopia::FormTraits<Left>::Implementation Implementation;
		typedef typename utopia::FormTraits<Left>::Scalar Scalar;

		static const int Order   = utopia::FormTraits<Left>::Order + FormTraits<Right>::Order;
		static const int Backend = utopia::FormTraits<Left>::Backend;
	};

	//variadic case
	template<class... Spaces>
	class Traits< ProductFunctionSpace<Spaces...> > {
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
		static const int Backend = utopia::FormTraits<First>::Backend;
	};
}

#endif //UTOPIA_PRODUCT_FUNCTION_SPACE_HPP
