#include "utopia_libmesh_ConvertTensor.hpp"

#include "utopia_Instance.hpp"

// All libmesh includes
#include "libmesh/petsc_vector.h"

namespace utopia {
    template <typename T>
    void ConvertTensor<UVector, libMesh::NumericVector<T>, 1>::apply(const UVector &in,
                                                                     libMesh::NumericVector<T> &out) {
#ifdef UTOPIA_ENABLE_PETSC
        auto ptr = dynamic_cast<libMesh::PetscVector<libMesh::Number> *>(&out);

        if (ptr) {
            Vec petsc_vec = ptr->vec();
            utopia::convert(in, petsc_vec);

        } else
#endif
        {
            assert(false && "IMPLEMENT ME");
            Utopia::Abort();
        }
    }

    template <typename T>
    void ConvertTensor<libMesh::NumericVector<T>, UVector, 1>::apply(const libMesh::NumericVector<T> &in,
                                                                     UVector &out) {
#ifdef UTOPIA_ENABLE_PETSC
        auto ptr = dynamic_cast<const libMesh::PetscVector<libMesh::Number> *>(&in);

        if (ptr) {
            Vec petsc_vec = ptr->vec();
            utopia::convert(petsc_vec, out);
        } else
#endif
        {
            assert(false && "IMPLEMENT ME");
            Utopia::Abort();
        }
    }

    template class ConvertTensor<libMesh::NumericVector<Traits<UVector>::Scalar>, UVector, 1>;
    template class ConvertTensor<UVector, libMesh::NumericVector<Traits<UVector>::Scalar>, 1>;
}  // namespace utopia
