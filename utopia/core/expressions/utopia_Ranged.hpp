#ifndef UTOPIA_UTOPIA_RANGED_HPP
#define UTOPIA_UTOPIA_RANGED_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Range.hpp"
#include "utopia_View.hpp"

namespace utopia {
    // template<class Derived, int Order>
    // class Ranged {};

    // template<class Derived>
    // class Ranged<Derived, 2> {
    // public:
      
    // };

    // template<class Derived>
    // class Ranged<Derived, 1> {
    // public:
      
    // };

    template<class Derived>
    inline View<Tensor<Derived, 2>> view(Tensor<Derived, 2> &t, const Range &row_range, const Range &col_range)
    {
        return View<Tensor<Derived, 2>>(t.derived(), row_range, col_range);
    }

    template<class Derived>
    inline View<Tensor<Derived, 1>> view(Tensor<Derived, 1> &t, const Range &range)
    {
        return View<Tensor<Derived, 1>>(t.derived(), range, Range(0));
    }
}

#endif //UTOPIA_UTOPIA_RANGED_HPP
