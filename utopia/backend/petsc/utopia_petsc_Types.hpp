
#ifndef UTOPIA_PETSC_TYPES_HPP
#define UTOPIA_PETSC_TYPES_HPP

#include "utopia_petsc_Traits.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_SparseMatrix.hpp"
#include "utopia_petsc_SerialSparseMatrix.hpp"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_SerialVector.hpp"
#include "utopia_petsc_Backend.hpp"

#include "utopia_Wrapper.hpp"

namespace utopia {

    UTOPIA_MAKE_TRAITS_DENSE(PetscMatrix, PetscTraits);
    UTOPIA_MAKE_TRAITS_SPARSE(PetscSparseMatrix, PetscTraits);
    UTOPIA_MAKE_TRAITS_SPARSE(PetscSerialSparseMatrix, PetscTraits);

    UTOPIA_MAKE_TRAITS(PetscVector, PetscTraits);
    UTOPIA_MAKE_TRAITS(PetscSerialVector, PetscTraits);

///////////////////////////////////////////////////////////////////////////////

    /*!
     * @brief Dense matrix representation of the petsc backend.
     */
    typedef utopia::Wrapper<PetscMatrix, 2>             DMatrixd;
    typedef utopia::Wrapper<PetscSparseMatrix, 2>       DSMatrixd;
    typedef utopia::Wrapper<PetscSerialSparseMatrix, 2> SSMatrixd;
    typedef utopia::Wrapper<PetscVector, 1>             DVectord;
    typedef utopia::Wrapper<PetscSerialVector, 1>       SVectord;

///////////////////////////////////////////////////////////////////////////////

    inline void disp(const Wrapper<PetscMatrix, 2> &w) {
        w.implementation().describe();
    }

    inline void disp(const Wrapper<PetscSparseMatrix, 2> &w) {
        w.implementation().describe();
    }

    inline void disp(const Wrapper<PetscSerialSparseMatrix, 2> &w) {
        w.implementation().describe();
    }

    inline void disp(const Wrapper<PetscVector, 1> &w) {
        w.implementation().describe();
    }

    inline void disp(const Wrapper<PetscSerialVector, 1> &w) {
        w.implementation().describe();
    }

///////////////////////////////////////////////////////////////////////////////

    /**
     * @brief      Unwrapps backend specific type of the tensor. \n
     *             Same routine applies to other utopia tensor types. 
     * @ingroup     interoperability
     *
     * @param      utopiaType  The wrapped/utopia type of tensor.
     */
    inline Mat &raw_type(Wrapper<PetscMatrix, 2> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline Mat &raw_type(Wrapper<PetscSparseMatrix, 2> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }


    inline Mat &raw_type(Wrapper<PetscSerialSparseMatrix, 2> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline Vec &raw_type(Wrapper<PetscVector, 1> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline Vec &raw_type(Wrapper<PetscSerialVector, 1> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

///////////////////////////////////////////////////////////////////////////////

    inline const Mat &raw_type(const Wrapper<PetscSparseMatrix, 2> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline const Mat &raw_type(const Wrapper<PetscSerialSparseMatrix, 2> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline const Mat &raw_type(const Wrapper<PetscMatrix, 2> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline const Vec &raw_type(const Wrapper<PetscVector, 1> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline const Vec &raw_type(const Wrapper<PetscSerialVector, 1> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }


    inline bool empty(const Wrapper<PetscVector, 1> &w)
    {
        return !w.implementation().initialized();
    }

    inline bool empty(const Wrapper<PetscSerialVector, 1> &w)
    {
        return !w.implementation().initialized();
    }

    inline void convert(const Vec &petsc_vec, Wrapper<PetscVector, 1> &utopia_vec)
    {
        Backend<PetscScalar, PETSC>::Instance().convert(petsc_vec, utopia_vec.implementation());
    }

    inline void convert(const Wrapper<PetscVector, 1> &utopia_vec, Vec &petsc_vec)
    {
        Backend<PetscScalar, PETSC>::Instance().convert(utopia_vec.implementation(), petsc_vec);
    }

    inline void convert(const Mat &petsc_mat, Wrapper<PetscMatrix, 2> &utopia_mat)
    {
        Backend<PetscScalar, PETSC>::Instance().convert(petsc_mat, utopia_mat.implementation());
    }

    inline void convert(const Wrapper<PetscSparseMatrix, 2> &utopia_mat, Mat &petsc_mat)
    {
        Backend<PetscScalar, PETSC>::Instance().convert(utopia_mat.implementation(), petsc_mat);
    }

    inline void convert(const Mat &petsc_mat, Wrapper<PetscSparseMatrix, 2> &utopia_mat)
    {
        Backend<PetscScalar, PETSC>::Instance().convert(petsc_mat, utopia_mat.implementation());
    }

    inline DSMatrixd sparse_mref(Mat &m) {
        DSMatrixd ret;
        ret.implementation().wrap(m);
        return ret;
    }

    inline DSMatrixd sparse_mref(const Mat &m) {
        DSMatrixd ret;
        ret.implementation().wrap(const_cast<Mat &>(m));
        return ret;
    }

    inline std::shared_ptr <DSMatrixd> sparse_mref_ptr(Mat &m) {
        DSMatrixd ret;
        Backend<PetscScalar, PETSC>::Instance().wrap(m, ret.implementation());
        return std::make_shared<DSMatrixd>(ret);
    }

    inline DVectord vref(Vec &vec)
    {
        DVectord ret;
        Backend<PetscScalar, PETSC>::Instance().wrap(vec, ret.implementation());
        return ret;
    }


///////////////////////////////////////////////////////////////////////////////
    
    UTOPIA_MAKE_PARALLEL_TRAITS(DMatrixd);
    UTOPIA_MAKE_PARALLEL_TRAITS(DVectord);
    UTOPIA_MAKE_PARALLEL_TRAITS(DSMatrixd);
}

#endif //UTOPIA_UTOPIA_PETSC_TYPES_HPP

