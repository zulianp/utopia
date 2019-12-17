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

    template<class VectorContainer>
    class VectorAccessor {
    public:
        using A = utopia::Accessor<VectorContainer>;

        template<typename T>
        static void set(VectorContainer &t, const std::size_t &i, const T &val1, const T &val2, const T &val3)
        {
            A::set(t, i, 0, val1);
            A::set(t, i, 1, val2);
            A::set(t, i, 2, val3);
        }

        template<typename T>
        static void set(VectorContainer &t, const std::size_t &i, const T &val1, const T &val2)
        {
            A::set(t, i, 0, val1);
            A::set(t, i, 1, val2);
        }
    };


}

#endif //UTOPIA_ACCESSOR_HPP