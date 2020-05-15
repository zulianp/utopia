#ifndef UTOPIA_FE_BASIC_OVERLOADS_HPP
#define UTOPIA_FE_BASIC_OVERLOADS_HPP

#include <type_traits>

namespace utopia {

    inline static double inner(const double left, const double right) { return left * right; }

    inline static float inner(const float left, const float right) { return left * right; }

    inline static int inner(const int left, const int right) { return left * right; }

    template <class T>
    struct remove_ref_and_const {
        typedef typename std::remove_const<typename std::remove_reference<T>::type>::type type;
    };
}  // namespace utopia

#endif  // UTOPIA_FE_BASIC_OVERLOADS_HPP
