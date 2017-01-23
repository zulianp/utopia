//
// Created by Patrick Zulian on 26/05/15.
//

#ifndef UTOPIA_UTOPIA_RANGED_HPP
#define UTOPIA_UTOPIA_RANGED_HPP

#include "utopia_Range.hpp"
#include "utopia_View.hpp"

namespace utopia {
    template<class Derived, int Order>
    class Ranged {};

    template<class Derived>
    class Ranged<Derived, 2> {
    public:
        View<Derived> range(const int rbegin, const int rend, const int cbegin, const int cend)
        {
            assert(rend - rbegin > 0);
            assert(rbegin >= 0);
            assert(rend > 0);

            assert(cend - cbegin > 0);
            assert(cbegin >= 0);
            assert(cend > 0);

            return View<Derived>(derived(), Range(rbegin, rend), Range(cbegin, cend));
        }

    private:
        DERIVED_CRT(Derived);
        CONST_DERIVED_CRT(Derived);
    };

    template<class Derived>
    class Ranged<Derived, 1> {
    public:
        View<Derived> range(const int begin, const int end)
        {
            assert(end - begin > 0);
            assert(begin >= 0);
            assert(end > 0);
            return View<Derived>(derived(), Range(begin, end), Range(0));
        }

    private:
        DERIVED_CRT(Derived);
        CONST_DERIVED_CRT(Derived);
    };
}

#endif //UTOPIA_UTOPIA_RANGED_HPP
