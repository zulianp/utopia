#ifndef UTOPIA_VIEW_FORWARD_DECLARATIONS_HPP
#define UTOPIA_VIEW_FORWARD_DECLARATIONS_HPP

#include <cassert>
//FIXME check if there is a device compilation
#define UTOPIA_DEVICE_ASSERT(expr) assert(expr)
// #define UTOPIA_DEVICE_ASSERT(...)

namespace utopia {
    using Size_t = std::size_t;
    static const Size_t DYNAMIC_SIZE = 0;

    template<class InnerExpr, class Op>
    class DeviceUnary;

    template<class InnerExpr>
    using DeviceNegate = utopia::DeviceUnary<InnerExpr, Minus>;

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

    template<class Left, class Right>
    class DeviceOuterProduct;

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

    template<class Expr>
    class DeviceTrace;

    template<class Expr>
    class DeviceDiag;

    template<class Expr>
    class DeviceEigenValues;


    template<class Expr>
    class DeviceEigenVectors;

    template<class Expr>
    class DeviceSingularValues;

    template<class T>
    class DeviceIdentity;

    template<typename T, Size_t... Args>
    class ArrayView;

    template<class View, int Order>
    class TensorView;

    template<class Expr, int Order, int TensorOrder>
    class DeviceNorm;

    template<class Expr, class Op, int Order>
    class DeviceReduce;

    template<class Expr, class Op, int Dim>
    class DeviceTensorReduce;

    template<class Left, class Right>
    class DeviceCrossProduct;

    template<class View>
    using VectorView = utopia::TensorView<View, 1>;

    template<class View>
    using MatrixView = utopia::TensorView<View, 2>;

    template<typename T, Size_t Size>
    using StaticVector = utopia::VectorView<ArrayView<T, Size>>;

    template<typename T>
    using StaticVector2 = utopia::StaticVector<T, 2>;

    template<typename T>
    using StaticVector3 = utopia::StaticVector<T, 3>;

    template<typename T>
    using StaticVector4 = utopia::StaticVector<T, 4>;

    template<typename T, Size_t Rows, Size_t Cols>
    using StaticMatrix = utopia::MatrixView<ArrayView<T, Rows, Cols>>;

    template<typename T>
    using StaticMatrix2x2 = utopia::StaticMatrix<T, 2, 2>;

    template<typename T>
    using StaticMatrix3x3 = utopia::StaticMatrix<T, 3, 3>;

    template<typename T>
    using StaticMatrix4x4 = utopia::StaticMatrix<T, 4, 4>;

}

#endif
