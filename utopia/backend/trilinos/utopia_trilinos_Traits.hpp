#ifndef UTOPIA_TRILINOSTRAITS_HPP
#define UTOPIA_TRILINOSTRAITS_HPP

#include "utopia_Traits.hpp"

#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_SparseMatrix.hpp"

#include "utopia_Tpetra_Vector.hpp"

#include "utopia_Base.hpp"

namespace utopia
{
class TrilinosTraits
    {
    public:
        typedef double                           Scalar;
        typedef int                              SizeType;
        typedef utopia::TpetraMatrix             Matrix;
        typedef utopia::TpetraSparseMatrix       SparseMatrix;
        //typedef utopia::TpetraSerialSparseMatrix SerialSparseMatrix; //local matrix

        typedef utopia::TpetraVector             Vector;
        //typedef utopia::TpetraSerialVector SerialVector;

        enum
            {
            Backend = TRILINOS
            };
    };
}

#endif //UTOPIA_TRILINOSTRAITS_HPP
