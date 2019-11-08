#ifndef UTOPIA_KOKKOS_TRAITS_HPP
#define UTOPIA_KOKKOS_TRAITS_HPP

#include "utopia_ViewForwardDeclarations.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    template<typename T, Size_t... Args>
    class Traits< ArrayView<T, Args...> > {
    public:
        using Scalar = T;
        using SizeType = Size_t;

        static const int Backend = SERIAL_HOMEMADE;
        static const int FILL_TYPE = FillType::DENSE;
    };

    template<class View>
    class Traits< VectorView<View> > : public Traits<View> {};

    template<class View>
    class Traits< MatrixView<View> > : public Traits<View> {};

    template<typename T, Size_t Size>
    using StaticVector = utopia::VectorView<ArrayView<T, Size>>; 

    template<typename T>
    using Vector2 = utopia::StaticVector<T, 2>;

    template<typename T>
    using Vector3 = utopia::StaticVector<T, 3>;

    template<typename T>
    using Vector4 = utopia::StaticVector<T, 4>;

    template<typename T, Size_t Rows, Size_t Cols>
    using StaticMatrix = utopia::MatrixView<ArrayView<T, Rows, Cols>>; 

    template<typename T>
    using Matrix2x2 = utopia::StaticMatrix<T, 2, 2>;

    template<typename T>
    using Matrix3x3 = utopia::StaticMatrix<T, 3, 3>;

    template<typename T>
    using Matrix4x4 = utopia::StaticMatrix<T, 4, 4>;


    //////////////////////////////////////// Type system extensions //////////


    template<class Left, class Right>
    class MultiplyResultType {};

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

        // typedef typename TensorQuery<Traits, Expr::Order>::Type Type;
    };
}

#endif //UTOPIA_KOKKOS_TRAITS_HPP
