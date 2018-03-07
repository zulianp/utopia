#ifndef UTOPIA_TRILINOSTRAITS_HPP
#define UTOPIA_TRILINOSTRAITS_HPP

#include "utopia_Traits.hpp"

//#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_SparseMatrix.hpp"

//#include "utopia_Tpetra_Vector.hpp"

//#include "utopia_Base.hpp"

namespace utopia
{
class TrilinosTraits
    {
    public:
        typedef double                   Scalar;
        typedef int                      SizeType;
        typedef TpetraMatrix             Matrix;
        typedef TpetraSparseMatrix       SparseMatrix;
        typedef TpetraVector             Vector;
        enum
            {
            Backend = TRILINOS
            };
    };
}

#endif //UTOPIA_TRILINOSTRAITS_HPP
