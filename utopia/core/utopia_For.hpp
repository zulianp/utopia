#ifndef UTOPIA_FOR_HPP
#define UTOPIA_FOR_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

	template<std::size_t Unroll = 32>
	class For {
	public:
		template<typename F>
		inline static void apply(
			const std::size_t &begin, 
			const std::size_t &end, 
			F f)
		{	
			const auto r = end - begin;
			auto i = begin;

			if(r < Unroll) {
				for(; i < end; ++i) {
					f(i);
				}

				return;
			}

			auto u_end = begin + (r/Unroll) * Unroll;

			for(; i < u_end; i += Unroll) {
				for(std::size_t k = 0; k < Unroll; ++k) {
					auto index = i + k;
					f(index);
				}
			}

			for(; i < end; ++i) {
				f(i);
			}
		}
	};
}

#endif //UTOPIA_FOR_HPP
