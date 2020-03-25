#ifndef UTOPIA_VIEW_TRAITS_HPP
#define UTOPIA_VIEW_TRAITS_HPP

#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template<typename T, Size_t... Args>
    class Traits< ArrayView<T, Args...> > {
    public:
        using Scalar = T;
        using SizeType = Size_t;
        using ValueType = T;

        static const int Backend = SERIAL_HOMEMADE;
        static const int FILL_TYPE = FillType::DENSE;
    };


    //////////////////////////////////////// Type system extensions //////////

    template<class Left, class Right>
    class MultiplyResultType {
    public:
        using Type = Right;
    };

    template<
        class LeftT, Size_t Rows, Size_t Cols,
        class RightT, Size_t RightSize
    >
    class MultiplyResultType<
        StaticMatrix<LeftT, Rows, Cols>,
        StaticVector<RightT, RightSize> > {
    public:
        using Type = utopia::StaticVector<RightT, Rows>;
    };

    template<class Traits, class View, int Order, int Sparsity>
    class TypeAndFill<Traits, Tensor<VectorView<View>, Order>, Sparsity> {
    public:
        using Type = utopia::VectorView<View>;
    };

    template<class Traits, class View, int Order, int Sparsity>
    class TypeAndFill<Traits, Tensor<MatrixView<View>, Order>, Sparsity> {
    public:
        using Type = utopia::MatrixView<View>;
    };

    template<class Traits, class Left, class Right, int Sparsity>
    class TypeAndFill<Traits, Multiply<Left, Right>, Sparsity> {
    public:
        using LeftType  = typename TypeAndFill<Traits, Left, Sparsity>::Type;
        using RightType = typename TypeAndFill<Traits, Right, Sparsity>::Type;
        using Type = typename MultiplyResultType<LeftType, RightType>::Type;
    };
}

#endif //UTOPIA_VIEW_TRAITS_HPP
