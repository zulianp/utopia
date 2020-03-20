#ifndef UTOPIA_ELEM_TRAITS_HPP
#define UTOPIA_ELEM_TRAITS_HPP

namespace utopia {
    template<class Elem>
    struct is_simplex {
        static const bool value = false;
    };
}

#endif //UTOPIA_ELEM_TRAITS_HPP
