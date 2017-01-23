
#ifndef UTOPIA_PETSC_TYPES_HPP
#define UTOPIA_PETSC_TYPES_HPP

#include "utopia_PETScTraits.hpp"
#include "utopia_PETScMatrix.hpp"
#include "utopia_PETScSparseMatrix.hpp"
#include "utopia_PETScSerialSparseMatrix.hpp"
#include "utopia_PETScVector.hpp"
#include "utopia_PETScSerialVector.hpp"
#include "utopia_PETScBackend.hpp"

#include "utopia_Wrapper.hpp"

namespace utopia {

    UTOPIA_MAKE_TRAITS_DENSE(PETScMatrix, PETScTraits);
    UTOPIA_MAKE_TRAITS_SPARSE(PETScSparseMatrix, PETScTraits);
    UTOPIA_MAKE_TRAITS_SPARSE(PETScSerialSparseMatrix, PETScTraits);

    UTOPIA_MAKE_TRAITS(PETScVector, PETScTraits);
    UTOPIA_MAKE_TRAITS(PETScSerialVector, PETScTraits);

///////////////////////////////////////////////////////////////////////////////

    /*!
     * @brief Dense matrix representation of the petsc backend.
     */
    typedef utopia::Wrapper<PETScMatrix, 2>             DMatrixd;
    typedef utopia::Wrapper<PETScSparseMatrix, 2>       DSMatrixd;
    typedef utopia::Wrapper<PETScSerialSparseMatrix, 2> SSMatrixd;
    typedef utopia::Wrapper<PETScVector, 1>             DVectord;
    typedef utopia::Wrapper<PETScSerialVector, 1>       SVectord;

///////////////////////////////////////////////////////////////////////////////

    inline void disp(const Wrapper<PETScMatrix, 2> &w) {
        w.implementation().describe();
    }

    inline void disp(const Wrapper<PETScSparseMatrix, 2> &w) {
        w.implementation().describe();
    }

    inline void disp(const Wrapper<PETScSerialSparseMatrix, 2> &w) {
        w.implementation().describe();
    }

    inline void disp(const Wrapper<PETScVector, 1> &w) {
        w.implementation().describe();
    }

    inline void disp(const Wrapper<PETScSerialVector, 1> &w) {
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
    inline Mat &raw_type(Wrapper<PETScMatrix, 2> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline Mat &raw_type(Wrapper<PETScSparseMatrix, 2> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }


    inline Mat &raw_type(Wrapper<PETScSerialSparseMatrix, 2> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline Vec &raw_type(Wrapper<PETScVector, 1> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline Vec &raw_type(Wrapper<PETScSerialVector, 1> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

///////////////////////////////////////////////////////////////////////////////

    inline const Mat &raw_type(const Wrapper<PETScSparseMatrix, 2> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline const Mat &raw_type(const Wrapper<PETScSerialSparseMatrix, 2> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline const Mat &raw_type(const Wrapper<PETScMatrix, 2> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline const Vec &raw_type(const Wrapper<PETScVector, 1> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }

    inline const Vec &raw_type(const Wrapper<PETScSerialVector, 1> &utopiaType)
    {
        return utopiaType.implementation().implementation();
    }


    inline bool empty(const Wrapper<PETScVector, 1> &w)
    {
        return !w.implementation().isInitialized();
    }

    inline bool empty(const Wrapper<PETScSerialVector, 1> &w)
    {
        return !w.implementation().isInitialized();
    }

    inline void convert(const Vec &petsc_vec, Wrapper<PETScVector, 1> &utopia_vec)
    {
        Backend<PetscScalar, PETSC>::Instance().convert(petsc_vec, utopia_vec.implementation());
    }

    inline void convert(const Mat &petsc_mat, Wrapper<PETScMatrix, 2> &utopia_mat)
    {
        Backend<PetscScalar, PETSC>::Instance().convert(petsc_mat, utopia_mat.implementation());
    }

    inline void convert(const Mat &petsc_mat, Wrapper<PETScSparseMatrix, 2> &utopia_mat)
    {
        Backend<PetscScalar, PETSC>::Instance().convert(petsc_mat, utopia_mat.implementation());
    }

    inline DSMatrixd sparse_mref(Mat &m) {
        DSMatrixd ret;
        Backend<PetscScalar, PETSC>::Instance().wrap(m, ret.implementation());
        return ret;
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

