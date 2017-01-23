//
// Created by Patrick Zulian on 26/05/15.
//

#ifndef UTOPIA_UTOPIA_STRUCTURED_HPP
#define UTOPIA_UTOPIA_STRUCTURED_HPP

#include "utopia_Size.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {
    template<class Derived>
    class Structured {
    public:
        inline Size size() const
        {
            using std::move;
            Size size;
            Evaluator<typename Derived::Implementation, Traits<Derived>::Backend> eval;
            eval.eval(derived(), size);
            return move(size);
        }

        inline void resize(const Size &size)
        {
            backend(this->derived()).resize(size, derived().implementation());
        }

        // template<int FOrder>
        // inline Derived &operator=(const Factory<Resize, FOrder> &f)
        // {
        //     resize(f.size());
        //     return derived();
        // }

    private:
        ALL_DERIVED_CRT(Derived)
    };

    template<class Derived, int Order>
    inline Size size(const Wrapper<Derived, Order> &wrapper)
    {
        return wrapper.size();
    }

    template<class Derived, int Order>
    inline Size local_size(const Wrapper<Derived, Order> &wrapper)
    {
        Size ret;
        backend(wrapper).local_size(wrapper.implementation(), ret);
        return ret;
    }

}

#endif //UTOPIA_UTOPIA_STRUCTURED_HPP