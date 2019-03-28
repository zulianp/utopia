
#ifndef UTOPIA_PETSC_TYPES_HPP
#define UTOPIA_PETSC_TYPES_HPP

#include "utopia_petsc_Traits.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_SparseMatrix.hpp"
#include "utopia_petsc_SerialSparseMatrix.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_SerialVector.hpp"
#include "utopia_petsc_Backend.hpp"
#include "utopia_Logger.hpp"

#include "utopia_Wrapper.hpp"

namespace utopia {
    UTOPIA_MAKE_TRAITS_DENSE(PetscMatrix, PetscTraits);
    UTOPIA_MAKE_TRAITS_SPARSE(PetscSparseMatrix, PetscTraits);
    UTOPIA_MAKE_TRAITS_SPARSE(PetscSerialSparseMatrix, PetscTraits);
    UTOPIA_MAKE_TRAITS_SPARSE(PetscCuSparseMatrix, PetscCudaTraits);

    UTOPIA_MAKE_TRAITS(PetscVector, PetscTraits);
    UTOPIA_MAKE_TRAITS(PetscSerialVector, PetscTraits);
    UTOPIA_MAKE_TRAITS(PetscCuVector, PetscCudaTraits);

    ///////////////////////////////////////////////////////////////////////////////

    /*!
     * @brief Dense matrix representation of the petsc backend.
     */
    using DMatrixd  = utopia::Wrapper<PetscMatrix, 2>;
    using DSMatrixd = utopia::Wrapper<PetscSparseMatrix, 2>;
    using SSMatrixd = utopia::Wrapper<PetscSerialSparseMatrix, 2>;
    using DVectord  = utopia::Wrapper<PetscVector, 1>;
    using SVectord  = utopia::Wrapper<PetscSerialVector, 1>;
    using CuSMatrixd = utopia::Wrapper<PetscCuSparseMatrix, 2>;
    using CuVectord  = utopia::Wrapper<PetscCuVector, 1>;

    /////////////////////////////////////////////////////////////////////////////////

    template<class Vector>
    inline void convert(const Vec &petsc_vec, Wrapper<Vector, 1> &utopia_vec)
    {
        Backend<PetscScalar, PETSC>::Instance().convert(petsc_vec, utopia_vec.implementation());
    }

    template<class Vector>
    inline void convert(const Wrapper<Vector, 1> &utopia_vec, Vec &petsc_vec)
    {
        Backend<PetscScalar, PETSC>::Instance().convert(utopia_vec.implementation(), petsc_vec);
    }

    template<class Matrix>
    inline void convert(const Mat &petsc_mat, Wrapper<Matrix, 2> &utopia_mat)
    {
        Backend<PetscScalar, PETSC>::Instance().convert(petsc_mat, utopia_mat.implementation());
    }

    template<class Matrix>
    inline void convert(const Wrapper<Matrix, 2> &utopia_mat, Mat &petsc_mat)
    {
        Backend<PetscScalar, PETSC>::Instance().convert(utopia_mat.implementation(), petsc_mat);
    }

    template<class Matrix>
    inline void wrap(Mat &m, Wrapper<Matrix, 2> &utopia_mat) {
        utopia_mat.implementation().wrap(m);
    }

    template<class Vector>
    inline void wrap(Vec &vec, Wrapper<Vector, 1> &utopia_vec)
    {
        Backend<PetscScalar, PETSC>::Instance().wrap(vec, utopia_vec.implementation());
    }

    inline int comm_size(const DVectord &t)
    {
        auto comm = t.implementation().communicator();
        int ret;
        MPI_Comm_size(comm, &ret);
        return ret;
    }

    inline int comm_rank(const DVectord &t)
    {
        auto comm = t.implementation().communicator();
        int ret;
        MPI_Comm_rank(comm, &ret);
        return ret;
    }

    inline int comm_size(const DSMatrixd &t)
    {
        auto comm = t.implementation().communicator();
        int ret;
        MPI_Comm_size(comm, &ret);
        return ret;
    }

    inline int comm_rank(const DSMatrixd &t)
    {
        auto comm = t.implementation().communicator();
        int ret;
        MPI_Comm_rank(comm, &ret);
        return ret;
    }

    inline void synchronize(DVectord &t)
    {
        t.implementation().update_ghosts();
    }


    UTOPIA_DEPRECATED_MSG("sparse_mref will be removed use wrap(Mat, DSMatrixd) instead")
    inline DSMatrixd sparse_mref(const Mat &m) {
        m_utopia_warning_once("sparse_mref will be removed use wrap(Mat, DSMatrixd) instead");
        DSMatrixd ret;
        wrap(const_cast<Mat &>(m), ret);
        return ret;
    }

    UTOPIA_DEPRECATED_MSG("sparse_mref_ptr will be removed use wrap(Mat, DSMatrixd) instead")
    inline std::shared_ptr<DSMatrixd> sparse_mref_ptr(Mat &m) {
        m_utopia_warning_once("sparse_mref_ptr will be removed use wrap(Mat, DSMatrixd) instead");
        auto ret = std::make_shared<DSMatrixd>();
        wrap(m, *ret);
        return ret;
    }

    ///////////////////////////////////////////////////////////////////////////////

    UTOPIA_MAKE_PARALLEL_TRAITS(DMatrixd);
    UTOPIA_MAKE_PARALLEL_TRAITS(DVectord);
    UTOPIA_MAKE_PARALLEL_TRAITS(DSMatrixd);
    UTOPIA_MAKE_PARALLEL_TRAITS(CuSMatrixd);
    UTOPIA_MAKE_PARALLEL_TRAITS(PetscCuVector);
}

#endif //UTOPIA_UTOPIA_PETSC_TYPES_HPP

