#ifndef UTOPIA_VIEW_FORWARD_DECLARATIONS_HPP
#define UTOPIA_VIEW_FORWARD_DECLARATIONS_HPP

//FIXME check if there is a device compilation
#define UTOPIA_DEVICE_ASSERT(expr) assert(expr)

namespace utopia {
    using Size_t = std::size_t;
    static const Size_t DYNAMIC_SIZE = 0;

    template<class InnerExpr, class Op>
    class DeviceUnary;

    template<class Left, class Right>
    class DeviceAssign;

    template<class Left, class Right, class Op>
    class DeviceBinary;

    template<class Left, class Right, class Op, int Order>
    class DeviceInPlace;

    template<class Left, class Right, int Order>
    class DeviceApproxEqual;

    template<class Left, class Right>
    class DeviceMultiply;

    template<class Left, class Right, int Order>
    class DeviceDot;

    template<class Derived>
    class DeviceExpression;

    template<class Derived>
    class DeviceTranspose;

    template<class T>
    class DeviceNumber;

    template<class T>
    class DeviceDeterminant;

    template<class Expr>
    class DeviceInverse;

    template<typename T, Size_t... Args>
    class ArrayView;

    template<class View, int Order>
    class TensorView;

    template<class View>
    using VectorView = utopia::TensorView<View, 1>;

    template<class View>
    using MatrixView = utopia::TensorView<View, 2>;

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

}

#endif
