#ifndef UTOPIA_FOR_HPP
#define UTOPIA_FOR_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Range.hpp"

namespace utopia {

    template<std::size_t Unroll = 1>
    class For {
    public:
        template<typename F>
        inline static void apply(const Range &r, F f)
        {
            apply(r.begin(), r.end(), f);
        }

        template<typename F>
        inline static void apply(
            const std::size_t &begin,
            const std::size_t &end,
            F f)
        {

            for(auto i = begin; i < end; ++i) {
                f(i);
            }

            // const auto r = end - begin;
            // auto i = begin;

            // if(r < Unroll) {
            // 	for(; i < end; ++i) {
            // 		f(i);
            // 	}

            // 	return;
            // }

            // auto u_end = begin + (r/Unroll) * Unroll;

            // for(; i < u_end; i += Unroll) {
            // 	for(std::size_t k = 0; k < Unroll; ++k) {
            // 		auto idx = i + k;
            // 		f(idx);
            // 	}
            // }

            // for(; i < end; ++i) {
            // 	f(i);
            // }
        }
    };

    template<int Backend>
    class ParallelFor : public For<1> {};
}

#endif //UTOPIA_FOR_HPP
