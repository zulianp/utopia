#ifndef UTOPIA_VIEW_FORWARD_DECLARATIONS_HPP
#define UTOPIA_VIEW_FORWARD_DECLARATIONS_HPP

namespace utopia {
    using Size_t = std::size_t;
    static const Size_t DYNAMIC_SIZE = 0;

    template<typename T, Size_t... Args>
    class ArrayView;

    template<class ArrayView_>
    class VectorView;

    template<class ArrayView_>
    class MatrixView;
}

#endif
