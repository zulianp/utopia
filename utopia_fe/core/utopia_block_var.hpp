#ifndef UTOPIA_PARAMETER_HPP
#define UTOPIA_PARAMETER_HPP 

#include "utopia.hpp"

#include <vector>
#include <unordered_map>
#include <utility>  

namespace utopia {

	template<class Tensor>
	class BlockVar : public Expression< BlockVar<Tensor> > {
	public:
		typedef typename utopia::Traits<Tensor>::Scalar Scalar;

		enum {
		    Order = utopia::Traits<Tensor>::Order
		};

		BlockVar(
			const Tensor &default_value,
			const std::vector< std::pair<int, Tensor> > &init_list)
		: default_value_(default_value)
		{
			for(auto &p : init_list) {
				values_[p.first] = p.second;
			}
		}

		const Tensor &get(const int id) const
		{
			auto it = values_.find(id);
			if(it == values_.end()) {
				return default_value_;
			}

			return it->second;
		}


		inline std::string getClass() const override
		{
			return "BlockVar";
		}

		Tensor default_value_;
		std::unordered_map<int, Tensor> values_;
	};

	template<class Tensor>
	inline BlockVar<Tensor> block_var(
		const Tensor &default_value,
		const std::vector< std::pair<int, Tensor> > &init_list)
	{
		return BlockVar<Tensor>(default_value, init_list);
	}

	template<class Expr>
	class Traits< BlockVar<Expr> > : public Traits<Expr> {
	public:
		enum {
			FILL_TYPE = utopia::FillType::DENSE
		};
	};
}


#endif //UTOPIA_PARAMETER_HPP
