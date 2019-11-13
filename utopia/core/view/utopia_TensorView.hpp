#ifndef UTOPIA_TENSOR_VIEW_HPP
#define UTOPIA_TENSOR_VIEW_HPP

#include "utopia_Tensor.hpp"

namespace utopia {
    template<class ArrayView, int Order>
    class TensorView {};

    template<class ArrayView, int Order_>
    class Traits< TensorView<ArrayView, Order_> > : public Traits<ArrayView> {
    public:
        static const int Order = Order_;
    };


    template<class ArrayView, int Order>
    inline void disp(const TensorView<ArrayView, Order> &t)
    {
        t.describe();
    }

    //hacks until I find a way to avoid this
    class DelegateArgs {};
    static const DelegateArgs args__;
}

#endif
