
#ifndef UTOPIA_TRILINOS_TYPES_HPP
#define UTOPIA_TRILINOS_TYPES_HPP

#include "utopia_trilinos_Traits.hpp"
#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_SparseMatrix.hpp"
#include "utopia_Tpetra_Vector.hpp"
#include "utopia_trilinos_Backend.hpp"

#include "utopia_Wrapper.hpp"


namespace utopia
{

UTOPIA_MAKE_TRAITS_DENSE(TpetraMatrix, TrilinosTraits);
UTOPIA_MAKE_TRAITS_SPARSE(TpetraSparseMatrix, TrilinosTraits);
//    UTOPIA_MAKE_TRAITS_SPARSE(TpetraSerialSparseMatrix, TpetraTraits);

UTOPIA_MAKE_TRAITS(TpetraVector, TrilinosTraits);
//UTOPIA_MAKE_TRAITS(TpetraSerialVector, TpetraTraits);


///////////////////////////////////////////////////////////////////////////////

/*!
 * @brief Dense matrix representation of the trilinos backend.
 */


typedef utopia::Wrapper<TpetraMatrix, 2>             DMatrixd;
typedef utopia::Wrapper<TpetraSparseMatrix, 2>       DSMatrixd;
//    typedef utopia::Wrapper<TpetraSerialSparseMatrix, 2> SSMatrixd;
typedef utopia::Wrapper<TpetraVector, 1>             DVectord;
//    typedef utopia::Wrapper<TpetraSerialVector, 1>       SVectord;

///////////////////////////////////////////////////////////////////////////////

/*    inline void disp(const Wrapper<TpetraMatrix, 2> &w) {
        w.implementation().describe();
    }

    inline void disp(const Wrapper<TpetraSparseMatrix, 2> &w) {
        w.implementation().describe();
    }

    inline void disp(const Wrapper<TpetraSerialSparseMatrix, 2> &w) {
        w.implementation().describe();
    }

    inline void disp(const Wrapper<TpetraVector, 1> &w) {
        w.implementation().describe();
    }

    inline void disp(const Wrapper<TpetraSerialVector, 1> &w) {
        w.implementation().describe()*/

///////////////////////////////////////////////////////////////////////////////
UTOPIA_MAKE_PARALLEL_TRAITS(DMatrixd);
UTOPIA_MAKE_PARALLEL_TRAITS(DVectord);
UTOPIA_MAKE_PARALLEL_TRAITS(DSMatrixd);
}

#endif //UTOPIA_UTOPIA_TRILINOS_TYPES_HPP

