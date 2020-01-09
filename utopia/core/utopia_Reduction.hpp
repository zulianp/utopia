#ifndef UTOPIA_REDUCTION_HPP
#define UTOPIA_REDUCTION_HPP

#include "utopia_For.hpp"

namespace utopia {

    class Reduction {
    public:
        template<typename F, class Accumulator>
        inline static void apply(const Range &r, F f, Accumulator &acc)
        {
            apply(r.begin(), r.end(), f, acc);
        }

        template<typename F, class Accumulator>
        inline static void apply(
            const std::size_t &begin,
            const std::size_t &end,
            F f,
            Accumulator &acc)
        {
            for(auto i = begin; i < end; ++i) {
                acc += f(i);
            }
        }

    };

    template<int Backend>
    class ParallelReduce : public Reduction {};
}

#endif //UTOPIA_REDUCTION_HPP
