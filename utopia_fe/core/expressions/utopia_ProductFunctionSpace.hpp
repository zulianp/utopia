#ifndef UTOPIA_PRODUCT_FUNCTION_SPACE_HPP
#define UTOPIA_PRODUCT_FUNCTION_SPACE_HPP 

#include "utopia_FormTraits.hpp"
#include "utopia_FunctionSpace.hpp"
#include <utility>
// #include <tuple>
#include <memory>

namespace utopia {

	template<class Space>
	class ProductFunctionSpace : public FunctionSpace<ProductFunctionSpace<Space> > {
	public:

		ProductFunctionSpace(std::initializer_list<std::shared_ptr<Space> > spaces)
		: spaces_(spaces.size())
		{ 
			std::copy(spaces.begin(), spaces.end(), spaces_.begin());
			init_subspace_ids();
		}

		//FIXME remove me
		void init_subspace_ids()
		{
			this->each([](const int i, Space &space) {
				if(space.subspace_id() < 0) {
					space.set_subspace_id(i);
				}
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

		template<class Fun>
		void each(Fun fun) const
		{
			for(std::size_t i = 0; i < spaces_.size(); ++i) {
				fun(i, subspace(i));
			}
		}

		inline void add_subspace(const std::shared_ptr<Space> &space)
		{
			//FIXME remove me
			if(space->subspace_id() < 0) {
				space->set_subspace_id(spaces_.size());
			}

			spaces_.push_back(space);
		}

		bool is_valid(const bool verbose)
		{
			auto n = spaces_.size();
			for(std::size_t i = 1; i < n; ++i) {
				if(spaces_[i-1]->subspace_id() + 1 != spaces_[i]->subspace_id()) {
					
					if(verbose) {
						std::cout << "subspace_" << (i-1) << " and subspace_" << (i) << "must have sequential ids" << std::endl;
						std::cout << (spaces_[i-1]->subspace_id()) << " + 1 != " << (spaces_[i]->subspace_id()) << std::endl;
					}

					return false;
				}
			}

			return true;
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
	class Traits<ProductFunctionSpace<Space>> : public Traits<Space> {
	public:
		// typedef typename utopia::FormTraits<Space>::Implementation Implementation;
		// typedef typename utopia::FormTraits<Space>::Scalar Scalar;
		// static const int Order   = utopia::FormTraits<Space>::Order;
		// static const int Backend = utopia::FormTraits<Space>::Backend;
	};

	template<class Space, class AssemblyContext>
	class FunctionalTraits<ProductFunctionSpace<Space>, AssemblyContext> {
	public:
		inline static int type(const ProductFunctionSpace<Space> &expr, const AssemblyContext &ctx)
		{
			int ret = 0;
			expr.each([&ret, &ctx](const int i, const Space &s) {
				ret = std::max(FunctionalTraits<Space, AssemblyContext>::type(s, ctx), ret);
			});

			return ret;
		}
		
		inline static int order(const ProductFunctionSpace<Space> &expr, const AssemblyContext &ctx)
		{
			int ret = 0;
			expr.each([&ret, &ctx](const int i, const Space &s) {
				ret = std::max(FunctionalTraits<Space, AssemblyContext>::order(s, ctx), ret);
			});

			return ret;
		}
	};
}

#endif //UTOPIA_PRODUCT_FUNCTION_SPACE_HPP
