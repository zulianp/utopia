
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

UTOPIA_MAKE_TRAITS(TpetraVector, TrilinosTraits);


///////////////////////////////////////////////////////////////////////////////

/*!
 * @brief Dense matrix representation of the trilinos backend.
 */


typedef Wrapper<TpetraMatrix, 2>             DMatrixd;
typedef Wrapper<TpetraSparseMatrix, 2>       DSMatrixd;
typedef Wrapper<TpetraVector, 1>             DVectord;

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

