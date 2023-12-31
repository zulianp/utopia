#ifndef UTOPIA_COMPARABLE_HPP
#define UTOPIA_COMPARABLE_HPP

#include "utopia_Traits.hpp"

namespace utopia {
    template <class T>
    class Comparable {
    public:
        virtual ~Comparable() = default;

        using Scalar = typename utopia::Traits<T>::Scalar;
        virtual bool equals(const T &other, const Scalar &tol = 0.0) const = 0;
    };
}  // namespace utopia

#endif  // UTOPIA_COMPARABLE_HPP
