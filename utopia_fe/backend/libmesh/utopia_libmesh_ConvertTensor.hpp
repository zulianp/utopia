#ifndef UTOPIA_LIBMESH_CONVERT_TENSOR_HPP
#define UTOPIA_LIBMESH_CONVERT_TENSOR_HPP

#include "utopia_ConvertTensor.hpp"

#include "utopia_libmesh_ForwardDeclarations.hpp"

#include "utopia_fe_base.hpp"

namespace utopia {

    template <typename T>
    class ConvertTensor<UVector, libMesh::NumericVector<T>, 1> {
    public:
        static void apply(const UVector &in, libMesh::NumericVector<T> &out);
    };

    template <typename T>
    void convert(const UVector &in, libMesh::NumericVector<T> &out) {
        ConvertTensor<UVector, libMesh::NumericVector<T>, 1>::apply(in, out);
    }

    template <typename T>
    class ConvertTensor<libMesh::NumericVector<T>, UVector, 1> {
    public:
        static void apply(const libMesh::NumericVector<T> &in, UVector &out);
    };

    template <typename T>
    void convert(const libMesh::NumericVector<T> &in, UVector &out) {
        ConvertTensor<libMesh::NumericVector<T>, UVector, 1>::apply(in, out);
    }

}  // namespace utopia

#endif  // UTOPIA_LIBMESH_CONVERT_TENSOR_HPP