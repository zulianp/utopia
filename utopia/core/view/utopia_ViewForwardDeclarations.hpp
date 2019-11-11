#ifndef UTOPIA_VIEW_FORWARD_DECLARATIONS_HPP
#define UTOPIA_VIEW_FORWARD_DECLARATIONS_HPP

namespace utopia {
    using Size_t = std::size_t;
    static const Size_t DYNAMIC_SIZE = 0;

    template<typename T, Size_t... Args>
    class ArrayView;

    template<class View, int Order>
    class TensorView;
    
    template<class View>
    using VectorView = utopia::TensorView<View, 1>;

    template<class View>
    using MatrixView = utopia::TensorView<View, 2>;
}

#endif
