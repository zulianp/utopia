
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

// UTOPIA_MAKE_TRAITS_DENSE(TpetraMatrix, TrilinosTraits);
UTOPIA_MAKE_TRAITS_SPARSE(TpetraMatrix, TrilinosTraits);

UTOPIA_MAKE_TRAITS(TpetraVector, TrilinosTraits);


///////////////////////////////////////////////////////////////////////////////

/*!
 * @brief Dense matrix representation of the trilinos backend.
 */


// typedef Wrapper<TpetraMatrix, 2>           TMatrixd;
typedef utopia::TpetraMatrix       TSMatrixd;
typedef Wrapper<TpetraVector, 1>             TVectord;

///////////////////////////////////////////////////////////////////////////////

    // inline void disp(const Wrapper<TpetraMatrix, 2> &w) {
    //     w.implementation().describe();
    // }

    // inline void disp(const Wrapper<TpetraSparseMatrix, 2> &w) {
    //     w.implementation().describe();
    // }

    // // inline void disp(const Wrapper<TpetraSerialSparseMatrix, 2> &w) {
    // //     w.implementation().describe();
    // // }

    // inline void disp(const TVectord &w) {
    //     w.implementation().describe();
    // }

    // inline int comm_size(const TVectord &t)
    // {
    //     auto comm = t.implementation().communicator();
    //     return comm->getSize();
    // }

    // inline int comm_rank(const TSMatrixd &t)
    // {
    //     auto comm = t.implementation().communicator();
    //     return comm->getRank();
    // }

    // inline int comm_size(const TSMatrixd &t)
    // {
    //     auto comm = t.implementation().communicator();
    //     return comm->getSize();
    // }

    // inline int comm_rank(const TVectord &t)
    // {
    //     auto comm = t.implementation().communicator();
    //     return comm->getRank();
    // }

    // inline void synchronize(TVectord &t)
    // {
    //     t.implementation().update_ghosts();
    // }
    // inline void disp(const Wrapper<TpetraSerialVector, 1> &w) {
    //     w.implementation().describe()
    // }

///////////////////////////////////////////////////////////////////////////////
UTOPIA_MAKE_PARALLEL_TRAITS(TMatrixd);
UTOPIA_MAKE_PARALLEL_TRAITS(TVectord);
UTOPIA_MAKE_PARALLEL_TRAITS(TSMatrixd);
}

#endif //UTOPIA_UTOPIA_TRILINOS_TYPES_HPP

