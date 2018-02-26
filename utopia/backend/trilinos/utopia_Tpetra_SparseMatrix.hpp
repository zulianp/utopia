#ifndef UTOPIA_TPETRASPARSEAMTRIX_H
#define UTOPIA_TPETRASPARSEMATRIX_H

#include "utopia_trilinos.hpp"

namespace utopia
{

class TpetraSparseMatrix : public TpetraMatrix
    {
    public:
        virtual ~TpetraSparseMatrix() {}
    };


#endif //UTOPIA_TPETRASPARSEMATRIX_H
