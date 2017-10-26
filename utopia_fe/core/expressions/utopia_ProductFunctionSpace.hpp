#ifndef UTOPIA_PRODUCT_FUNCTION_SPACE_HPP
#define UTOPIA_PRODUCT_FUNCTION_SPACE_HPP 

#include "utopia_FormTraits.hpp"
#include "utopia_FunctionSpace.hpp"
#include <utility>
// #include <tuple>
#include <memory>

namespace utopia {

	// template<class Space, int Begin, int End>
	// class SubspaceIterator {
	// public:
	// 	template<class Fun>
	// 	void visit(Fun fun)
	// 	{
	// 		fun(Begin, space.template get<Begin>() );
	// 		SubspaceIterator<Space, Begin + 1, End> next(space);
	// 		next.visit(fun);
	// 	}

	// 	SubspaceIterator(Space &space) : space(space) {}
	// 	Space &space;
	// };

	// template<class Space, int Begin>
	// class SubspaceIterator<Space, Begin, Begin + 1> {
	// public:
	// 	template<class Fun>
	// 	void visit(Fun fun)
	// 	{
	// 		fun(Begin, space.template get<Begin>());
	// 	}

	// 	SubspaceIterator(Space &space) : space(space) {}
	// 	Space &space;
	// };

	// template<class... Spaces>
	// class StaticProductFunctionSpace : public FunctionSpace<StaticProductFunctionSpace<Spaces...> > {
	// public:

	// 	static const int NumberOfSubSpaces = std::tuple_size< std::tuple<Spaces...> >::value;

	// 	StaticProductFunctionSpace(const Spaces &...spaces)
	// 	: spaces_(spaces...)
	// 	{ }

	// 	template<int Index>
	// 	inline auto get() -> decltype(std::get<Index>(std::tuple<Spaces...>())) &
	// 	{
	// 		return std::get<Index>(this->spaces_);
	// 	}

	// 	inline constexpr static std::size_t size() 
	// 	{
	// 		return NumberOfSubSpaces;
	// 	}

	// 	template<class Fun>
	// 	void each(Fun fun)
	// 	{
	// 		SubspaceIterator<StaticProductFunctionSpace, 0, NumberOfSubSpaces> iter(*this);
	// 		iter.visit(fun);
	// 	}

	// private:
	// 	std::tuple<Spaces...> spaces_;
	// };



	// // template<class LeftSpace, class RightSpace>
	// // inline StaticProductFunctionSpace<LeftSpace, RightSpace> operator * (const FunctionSpace<LeftSpace> &left, const FunctionSpace<RightSpace> &right)
	// // {
	// // 	return StaticProductFunctionSpace<LeftSpace, RightSpace>(left.derived(), right.derived());
	// // }

	// template<class... Spaces>
	// inline StaticProductFunctionSpace<FunctionSpace<Spaces>...> static_kron_prod(const FunctionSpace<Spaces> &... spaces)
	// {
	// 	return StaticProductFunctionSpace<FunctionSpace<Spaces>...>(spaces...);
	// }

	// //base case
	// template<class Left, class Right>
	// class Traits<StaticProductFunctionSpace<Left, Right>> {
	// public:
	// 	typedef typename utopia::FormTraits<Left>::Implementation Implementation;
	// 	typedef typename utopia::FormTraits<Left>::Scalar Scalar;

	// 	static const int Order   = utopia::FormTraits<Left>::Order + FormTraits<Right>::Order;
	// 	static const int Backend = utopia::FormTraits<Left>::Backend;
	// };

	// //variadic case
	// template<class... Spaces>
	// class Traits< StaticProductFunctionSpace<Spaces...> > {
	// public:
	// 	template<class Space>
	// 	inline static constexpr int order()
	// 	{
	// 		return Traits<Space>::Order;
	// 	}

	// 	template<class First, class...Rest>
	// 	inline static constexpr int order()
	// 	{
	// 		return Traits<First>::Order + order<Rest...>();
	// 	}

	// 	template<class First, class...Rest>
	// 	class GetFirst {
	// 	public:
	// 		typedef First Type;
	// 	};

	// 	typedef typename GetFirst<Spaces...>::Type First;
	// 	typedef typename utopia::FormTraits<First>::Implementation Implementation;
	// 	typedef typename utopia::FormTraits<First>::Scalar Scalar;

	// 	const static int Order = order<Spaces...>();
	// 	static const int Backend = utopia::FormTraits<First>::Backend;
	// };


	template<class Space>
	class ProductFunctionSpace : public FunctionSpace<ProductFunctionSpace<Space> > {
	public:

		ProductFunctionSpace(std::initializer_list<std::shared_ptr<Space> > spaces)
		: spaces_(spaces.size())
		{ 
			std::copy(spaces.begin(), spaces.end(), spaces_.begin());
			init_subspace_ids();
		}

		void init_subspace_ids()
		{
			this->each([](const int i, Space &space) {
				space.set_subspace_id(i);
			});
		}

		inline Space &subspace(const int index)
		{
			return *spaces_[index];
		}

		inline Space &operator[](const int index)
		{
			return *spaces_[index];
		}


		inline const Space &subspace(const int index) const
		{
			return *spaces_[index];
		}

		inline const Space &operator[](const int index) const
		{
			return *spaces_[index];
		}

		inline std::shared_ptr<Space> &subspace_ptr(const int index)
		{
			return spaces_[index];
		}

		inline std::size_t n_subspaces() 
		{
			return spaces_.size();
		}

		template<class Fun>
		void each(Fun fun)
		{
			for(std::size_t i = 0; i < spaces_.size(); ++i) {
				fun(i, subspace(i));
			}
		}

		inline void add_subspace(const std::shared_ptr<Space> &space)
		{
			space->set_subspace_id(spaces_.size());
			spaces_.push_back(space);
		}

	private:
		std::vector< std::shared_ptr<Space> > spaces_;
	};

	template<class Space>
	inline ProductFunctionSpace<Space> operator * (const FunctionSpace<Space> &left, const FunctionSpace<Space> &right)
	{
		return ProductFunctionSpace<Space>({ std::make_shared<Space>(left.derived()), std::make_shared<Space>(right.derived())});
	}

	template<class Space>
	inline ProductFunctionSpace<Space> operator * (const FunctionSpace<ProductFunctionSpace<Space> > &left, const FunctionSpace<Space> &right)
	{
		 ProductFunctionSpace<Space> ret(left.derived());
		 ret.add_subspace(std::make_shared<Space>(right.derived()));
		 return ret;
	}

	//base case
	template<class Space>
	class Traits<ProductFunctionSpace<Space>> {
	public:
		typedef typename utopia::FormTraits<Space>::Implementation Implementation;
		typedef typename utopia::FormTraits<Space>::Scalar Scalar;
		static const int Order   = utopia::FormTraits<Space>::Order;
		static const int Backend = utopia::FormTraits<Space>::Backend;
	};
}

#endif //UTOPIA_PRODUCT_FUNCTION_SPACE_HPP
