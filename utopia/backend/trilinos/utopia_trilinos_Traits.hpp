#ifndef UTOPIA_TRILINOSTRAITS_HPP
#define UTOPIA_TRILINOSTRAITS_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_trilinos_ForwardDeclaration.hpp"

namespace utopia {
class TrilinosTraits
    {
    public:
        typedef TpetraMatrix::SC         Scalar;
        typedef TpetraMatrix::GO         SizeType;
        typedef TpetraMatrix             Matrix;
        // typedef TpetraSparseMatrix       SparseMatrix;
        typedef TpetraVector             Vector;
        enum {
            Backend = TRILINOS
        };
    };
}

#endif //UTOPIA_TRILINOSTRAITS_HPP
