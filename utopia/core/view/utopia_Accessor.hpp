#ifndef UTOPIA_ACCESSOR_HPP
#define UTOPIA_ACCESSOR_HPP

#include "utopia_Base.hpp"
#include <array>

namespace utopia {

    template<class Container>
    class Accessor {};

    template<typename T, std::size_t N>
    class Accessor<std::array<T, N>> {
    public:
        static const T &get(const std::array<T, N> &t, const std::size_t &i)
        {
            return t[i];
        }

        static void set(std::array<T, N> &t, const std::size_t &i, const T &val)
        {
            t[i] = val;
        }
    };

}

#endif //UTOPIA_ACCESSOR_HPP