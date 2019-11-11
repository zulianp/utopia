#ifndef UTOPIA_TENSOR_VIEW_HPP
#define UTOPIA_TENSOR_VIEW_HPP

#include "utopia_Tensor.hpp"

namespace utopia {
    template<class ArrayView, int Order>
    class TensorView /*: public Tensor<TensorView<ArrayView, Order>, Order>)*/ {};

    template<class ArrayView, int Order_>
    class Traits< TensorView<ArrayView, Order_> > : public Traits<ArrayView> {
    public:
        static const int Order = Order_;
    };
}

#endif
