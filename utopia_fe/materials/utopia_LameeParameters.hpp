#ifndef UTOPIA_LAMEE_PARAMETERS_HPP
#define UTOPIA_LAMEE_PARAMETERS_HPP 

#include "utopia_block_var.hpp"
#include <map>
#include <vector>
#include <utility>

namespace utopia {
	class LameeParameters {
	public:
		LameeParameters(const double default_mu = 1, const double default_lambda = 1)
		: default_mu(default_mu), default_lambda(default_lambda)
		{}

		void set_mu(const int id, const double mu)
		{
			mu_[id] = mu;
		}
		
		void set_lambda(const int id, const double lambda)
		{
			lambda_[id] = lambda;
		}
		
		double mu(const int id) const
		{
			auto it = mu_.find(id);
			if(it == mu_.end()) {
				return default_mu;
			}
			
			return it->second;
		}
		
		double lambda(const int id) const
		{
			auto it = lambda_.find(id);
			if(it == lambda_.end()) {
				return default_lambda;
			}
			
			return it->second;
		}

		static BlockVar<double> make_var_from_map(const double default_val, const std::map<int, double> &pairs)
		{
			std::vector<std::pair<int, double> > vec;
			vec.reserve(pairs.size());

			for(const auto &p : pairs) {
				vec.push_back(p);
			}

			return BlockVar<double>(default_val, vec);
		}

		inline BlockVar<double> var_lambda() const
		{
			return make_var_from_map(default_lambda, lambda_);
		}

		inline BlockVar<double> var_mu() const
		{
			return make_var_from_map(default_mu, mu_);
		}
		
		double default_mu, default_lambda;
		std::map<int, double> mu_;
		std::map<int, double> lambda_;
	};
}

#endif //UTOPIA_LAMEE_PARAMETERS_HPP
