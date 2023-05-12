#ifndef UTOPIA_FOR_HPP
#define UTOPIA_FOR_HPP

#include "utopia_Base.hpp"
#include "utopia_Range.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template <std::size_t Unroll = 1>
    class For {
    public:
        template <typename F>
        inline static void apply(const Range &r, F f) {
            apply(r.begin(), r.end(), f);
        }

        template <typename F>
        inline static void apply(const std::size_t &begin, const std::size_t &end, F f) {
            for (auto i = begin; i < end; ++i) {
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

    template <int Backend>
    class ParallelFor {
    public:
        template <typename F>
        inline static void apply(const Range &r, F f) {
            apply(r.begin(), r.end(), f);
        }

        template <class T, typename F>
        inline static void apply(const RangeDevice<T> &r, F f) {
            apply(r.begin(), r.end(), f);
        }

        template <typename F>
        inline static void apply(const std::size_t &begin, const std::size_t &end, F f) {
            For<1>::apply(begin, end, f);
        }
    };
}  // namespace utopia

#endif  // UTOPIA_FOR_HPP
